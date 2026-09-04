"""Pairwise BioPython aligner wrapped in the mmalignments Element model."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import Bio  # type: ignore[import]
from Bio.Align import Alignment, PairwiseAligner  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

from mmalignments.models.artifacts import (
    ArtifactSet,
)
from mmalignments.models.elements import (
    CallSpec,
    Element,
    FastqSource,
    FileSource,
    element,
)
from mmalignments.models.externals import Runnable
from mmalignments.models.overlay import (
    CfgType,
    ElementTag,
    Out,
    OutputSpec,
    OutType,
    Params,
    ParType,
    TagType,
)
from mmalignments.models.preprocess import (  # type: ignore[import]
    EnsemblAPI,
    ReferenceSpec,
)
from mmalignments.models.tags import (
    Method,
    Omics,
    Stage,
    State,
)
from mmalignments.services.io import (
    read_frame,
    read_sequence_file,
    write_frames,
)

logger = logging.getLogger(__name__)


def _overlay_signature_determinants(params: ParType | None) -> tuple[str, ...]:
    if params is None:
        return ()
    merged = params.to_dict()
    tokens: list[str] = []
    for key in sorted(merged.keys()):
        value = merged[key]
        if value is None or value is False or value == []:
            continue
        tokens.append(f"{key}={value}")
    return tuple(tokens)


@dataclass(frozen=True)
class Variant:
    """Simple variant representation generated from a pairwise alignment."""

    chrom: str
    pos: int
    end: int
    ref: str
    alt: str
    variant_type: str
    description: str

    @property
    def is_ambiguous(self) -> bool:
        """Whether the alternate allele contains an ambiguous base."""
        return "N" in self.alt.upper()


@dataclass
class AlignmentResult:
    alignment: Alignment | None
    status: str
    reference_coverage: float | None = None
    query_coverage: float | None = None

class BioAligner:
    """Wrapper around Bio.Align.PairwiseAligner.

    The class offers two element-producing methods:
    - align: writes an alignment table
    - align2variant: writes a one-variant-per-row table and an alignment table
    """

    def __init__(
        self,
        *,
        mode: str = "global",
        match_score: float = 2.0,
        mismatch_score: float = -2.0,
        gap_open_score: float = -5.0,
        gap_extend_score: float = -1.0,
        free_query_ends: bool = True,
        min_reference_coverage: float = 0.90,
        min_query_coverage: float = 0,
    ) -> None:
        self.mode = mode
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_open_score = gap_open_score
        self.gap_extend_score = gap_extend_score
        self.free_query_ends = free_query_ends
        self.min_reference_coverage = min_reference_coverage
        self.min_query_coverage = min_query_coverage
        self.aligner = self.create_aligner()
        self.version_name = f"bioaligner_{Bio.__version__}"

    def create_aligner(self) -> PairwiseAligner:
        aligner = PairwiseAligner()
        aligner.mode = self.mode
        aligner.match_score = self.match_score
        aligner.mismatch_score = self.mismatch_score
        aligner.open_gap_score = self.gap_open_score
        aligner.extend_gap_score = self.gap_extend_score
        if self.free_query_ends:
            aligner.query_left_open_gap_score = 0
            aligner.query_left_extend_gap_score = 0
            aligner.query_right_open_gap_score = 0
            aligner.query_right_extend_gap_score = 0
        return aligner

    def _align(
        self,
        reference_sequence: str,
        query_sequence: str,
    ) -> AlignmentResult:
        alignments = self.aligner.align(reference_sequence, query_sequence)

        if not alignments:
            return AlignmentResult(alignment=None, status="failed", reference_coverage=None, query_coverage=None)

        alignment = alignments[0]
        reference_coverage = self._reference_coverage(alignment, reference_sequence)
        query_coverage = self._query_coverage(alignment, query_sequence)
        if reference_coverage < self.min_reference_coverage:
            return AlignmentResult(alignment=None, status="low_reference_coverage", reference_coverage=reference_coverage, query_coverage=query_coverage)

        if query_coverage < self.min_query_coverage:
            return AlignmentResult(alignment=None, status="low_query_coverage", reference_coverage=reference_coverage, query_coverage=query_coverage)

        return AlignmentResult(alignment=alignment, status="success", reference_coverage=reference_coverage, query_coverage=query_coverage)

    def _query_coverage(
        self,
        alignment: Alignment,
        query_sequence: str,
    ) -> float:
        aligned_query = alignment[1]
        aligned_bases = sum(base != "-" for base in aligned_query)
        return aligned_bases / len(query_sequence)

    def _reference_coverage(
        self,
        alignment: Alignment,
        reference_sequence: str,
    ) -> float:
        aligned_reference = alignment[0]
        aligned_bases = sum(base != "-" for base in aligned_reference)
        return aligned_bases / len(reference_sequence)

    @property
    def default_path(self) -> Path:
        return Path("results/alignments") / self.version_name

    @property
    def default_align_params(self) -> Params:
        return Params(
            species="species",
            assembly="assembly",
            chromosome="chromosome",
            start="start",
            stop="stop",
            strand="strand",
            sleep=0.1,
            query_id="R1 Name",
            query_sequence="R1",
            target_id="target",
            target_sequence="amplicon",
        )

    @property
    def determinants(self) -> tuple[str, ...]:
        return (
            f"mode={self.mode}",
            f"match={self.match_score}",
            f"mismatch={self.mismatch_score}",
            f"gap_open={self.gap_open_score}",
            f"gap_extend={self.gap_extend_score}",
        )

    @element
    def align(
        self,
        sequences: Element | FileSource | FastqSource,
        reference: (
            Element | FileSource | ReferenceSpec | Mapping[str, ReferenceSpec]
        ),  # noqa: E501
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """
        Returns an Element that aligns queries to their reference sequence and
        saves the alignments in a DataFrame.

        Parameters
        ----------
        sequences : Element | FileSource | FastqSource
            Input sequences to align. Can be an Element, a FileSource, or a
            FastqSource.
        reference : Element | FileSource | ReferenceSpec | Mapping[str,
            ReferenceSpec]
            Reference sequences to align against. Can be an Element, a
            FileSource, or a Mapping[str, ReferenceSpec].
                tag : TagType | None
            Optional tag override.
        out : OutType | None
            Output specification.
        par : ParType | None
            Additional parameter overrides.
        cfg : CfgType |None
            Subprocess configuration.

        Returns
        -------
        Element
            _description_

        Raises
        ------
        ValueError
            _description_
        """

        if isinstance(sequences, (Element, FileSource)):
            sequence_paths = (sequences.file,)
        elif isinstance(sequences, FastqSource):
            sequence_paths = sequences.primary.files
        else:
            raise ValueError(
                "Unsupported sequences type. Use Element, FileSource, or FastqSource."  # noqa: E501
            )
        tag = sequences.tag.bump(
            stage=Stage.ANALYSIS,
            method=Method.BIO,
            state=State.PROCESSED,
            omics=Omics.DNA,
            flag="bioalign",
        ).resolve(tag)
        out = OutputSpec.from_tag(tag, Out(folder=self.default_path, ext="tsv"))
        artifacts = ArtifactSet.from_outspec(out)
        ref_source = reference
        if isinstance(reference, (Element, FileSource)):
            ref_source = reference.file

        par = self.default_align_params.resolve(par)
        runner = self.align_and_write(
            sequence_paths=sequence_paths,
            reference_source=ref_source,
            out_paths=out.files,
            par=par,
        )

        determinants = self.determinants + _overlay_signature_determinants(par)
        key = Element.generate_key(tag, "BioAligner", "align")
        pres = (
            (sequences,)
            if not isinstance(reference, (Element, FileSource))
            else (sequences, reference)
        )
        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            pres=pres,
        )

    @element
    def align2variant(
        self,
        sequences: Element | FileSource | FastqSource,
        reference: (
            Element | FileSource | ReferenceSpec | Mapping[str, ReferenceSpec]
        ),  # noqa: E501
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """
        Returns an Element that aligns queries to their reference sequence and
        saves the alignments in a DataFrame.
        """

        if isinstance(sequences, (Element, FileSource)):
            sequence_paths = (sequences.file,)
        elif isinstance(sequences, FastqSource):
            sequence_paths = sequences.primary.files
        else:
            raise ValueError(
                "Unsupported sequences type. Use Element, FileSource, or FastqSource."  # noqa: E501
            )
        tag = sequences.tag.bump(
            stage=Stage.ANALYSIS,
            method=Method.BIO,
            state=State.PROCESSED,
            omics=Omics.DNA,
            flag="alignvariant",
        ).resolve(tag)
        out = OutputSpec.from_tag(tag, Out(folder=self.default_path, ext="tsv"))
        alignment_file = out.file.with_suffix(".alignments.tsv")
        artifacts = ArtifactSet.from_outspec(out).with_extra(
            "alignments", alignment_file
        )
        ref_source = reference
        if isinstance(reference, (Element, FileSource)):
            ref_source = reference.file

        runner = self.align_and_annotate_variant(
            sequence_paths=sequence_paths,
            reference_source=ref_source,
            variants_out=out.files,
            alignments_out=[alignment_file],
            par=par,
        )

        determinants = self.determinants + _overlay_signature_determinants(par)
        key = Element.generate_key(tag, "BioAligner", "align")
        pres = (
            (sequences,)
            if not isinstance(reference, (Element, FileSource))
            else (sequences, reference)
        )
        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            pres=pres,
        )

    @element
    def fetchref(
        self,
        reference: Element | FileSource | ReferenceSpec | Mapping[str, ReferenceSpec],
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """
        Fetches the reference sequences if they are not in the frame.
        """
        par = self.default_align_params.resolve(par)
        if isinstance(reference, (Element, FileSource)):
            ref_source = reference.file
            pres = (reference,)
            tag = reference.tag.bump(
                stage=Stage.ANALYSIS,
                method=Method.BIO,
                state=State.PROCESSED,
                flag="fetchref",
            ).resolve(tag)
        else:
            ref_source = reference
            pres = ()
            tag = ElementTag(
                root="references",
                level=0,
                stage=Stage.ANALYSIS,
                method=Method.BIO,
                state=State.PROCESSED,
                omics=Omics.DNA,
                flag="fetchref",
            )

        out = OutputSpec.from_tag(tag, Out(folder=Path("cache/references"), ext="tsv"))
        artifacts = ArtifactSet.from_outspec(out)
        runner = self.fetch_references(
            reference_source=ref_source,
            outfiles=out.files,
            par=par,
        )
        determinants = self.determinants + _overlay_signature_determinants(par)
        key = Element.generate_key(tag, "BioAligner", "align")
        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            pres=pres,
        )

    def fetch_references(
        self,
        reference_source: Path | ReferenceSpec | Mapping[str, ReferenceSpec],
        outfiles: list[Path],
        par: Params,
    ) -> Runnable:
        """Resolve references from file, mapping, literal sequence, or spec."""

        def run():
            frame = self.resolve_references(
                reference_source=reference_source,
                par=par,
            )
            write_frames(frame, outfiles)

        callspec = CallSpec(
            path=("BioAligner", "fetch_references"),
            args=(reference_source, outfiles),
            kwargs={"par": par},
        )
        return Runnable(run, display=callspec.render())

    def resolve_sequences(
        self, sequence_paths: tuple[Path, ...], par: Params
    ) -> DataFrame:
        """Resolve sequences from file, mapping, literal sequence, or spec."""
        if len(sequence_paths) == 1:
            return read_sequence_file(sequence_paths[0])
        elif len(sequence_paths) == 2:  # paired
            df1 = read_sequence_file(sequence_paths[0])
            df2 = read_sequence_file(sequence_paths[1])
            if len(df1) != len(df2):
                raise ValueError(
                    f"Paired sequence files have different lengths: {len(df1)} vs {len(df2)}"  # noqa: E501
                )
            df1["R2"] = df2.iloc[:, 0].values
            return df1
        else:
            raise ValueError(
                "Unsupported sequence source type. Use Element/FileSource/Path/str."  # noqa: E501
            )

    def annotate_alignment(
        self,
        sequences: DataFrame,
        references: DataFrame,
        alignments: Mapping[str, AlignmentResult],
        *,
        par: ParType | None = None,
    ) -> DataFrame:
        """
        Annotates calculated alignments. The resulting alignments are saved in a
        DataFrame.

        This columns in the returned frame are:

        `sequence_id` `sequence_column` `reference_column` score query_aligned
        reference_aligned

        Parameters
        ----------
        sequences : DataFrame
            _description_
        reference : Mapping[str, str]
            _description_
        par : ParType | None, optional
            _description_, by default None

        Returns
        -------
        DataFrame
            _description_
        """
        par = self.default_align_params.resolve(par)
        score_col = "score"
        reference_coverage_col = "reference_coverage"
        query_coverage_col = "query_coverage"
        qaligned_col = "query_aligned"
        taligned_col = "target_aligned"
        status_col = "alignment_status"
        # target_seq_col = par.get("target_sequence", "sequence")
        target_id_col = par.get("target_id", "target_id")
        query_id_col = par.get("query_id", "query_id")
        chromosome_col = par.get("chromosome", "chromosome")
        start_col = par.get("start", "start")
        stop_col = par.get("stop", "stop")
        strand_col = par.get("strand", "strand")
        species_col = par.get("species", "species")
        assembly_col = par.get("assembly", "assembly")
        ref_cols = [
            chromosome_col,
            start_col,
            stop_col,
            strand_col,
            species_col,
            assembly_col,
        ]
        additional_columns = {
            key: [] for key in [score_col, qaligned_col, taligned_col, status_col, reference_coverage_col, query_coverage_col] + ref_cols
        }
        for _, row in sequences.iterrows():
            target_id = row[target_id_col]
            reference_row = references.loc[references[target_id_col] == target_id].iloc[
                0
            ]
            query_id = row[query_id_col]
            alignment_result = alignments[query_id]
            target_aligned_seq, query_aligned_seq = "", ""
            score = float("-inf")
            if alignment_result.alignment is not None:
                target_aligned_seq, query_aligned_seq = alignment_to_strings(alignment_result.alignment)
                score = alignment_result.alignment.score
            additional_columns[score_col].append(score)
            additional_columns[qaligned_col].append(query_aligned_seq)
            additional_columns[taligned_col].append(target_aligned_seq)
            additional_columns[status_col].append(alignment_result.status)
            additional_columns[reference_coverage_col].append(alignment_result.reference_coverage)
            additional_columns[query_coverage_col].append(alignment_result.query_coverage)
            for col in ref_cols:
                additional_columns[col].append(reference_row[col])

        df = sequences.copy()
        for col in additional_columns:
            df[col] = additional_columns[col]

        return df

    def align_and_write(
        self,
        sequence_paths: tuple[Path, ...],
        reference_source: Path,
        out_paths: list[Path],
        par: ParType | None = None,
    ) -> Runnable:
        def run():
            params = self.default_align_params.resolve(par)
            references = read_frame(reference_source)
            sequences = self.resolve_sequences(sequence_paths, par=params)
            sequences = self.annotate_ns(sequences, par=params)
            alignments = self.align_sequences(sequences, references, par=params)
            annotated = self.annotate_alignment(
                sequences, references, alignments, par=params
            )
            write_frames(annotated, out_paths)
            return annotated

        callspec = CallSpec(
            path=("BioAligner", "align_and_write"),
            args=(sequence_paths, reference_source, out_paths),
            kwargs={"par": par},
        )
        return Runnable(run, display=callspec.render())

    def align_and_annotate_variant(
        self,
        sequence_paths: tuple[Path, ...],
        reference_source: Path,
        variants_out: list[Path],
        alignments_out: list[Path],
        *,
        par: ParType | None = None,
    ) -> Runnable:
        def __call() -> DataFrame:
            # start = time.perf_counter()
            params = self.default_align_params.resolve(par)
            # elapsed = time.perf_counter() - start
            # print(f"Elapsed time for resolving parameters: {elapsed:.4f} seconds")
            references = read_frame(reference_source)
            # print(f"Elapsed time for resolving references: {elapsed:.4f} seconds")
            sequences = self.resolve_sequences(sequence_paths, par=params)
            sequences = self.annotate_ns(sequences, par=params)
            # print(sequences.head())
            # elapsed = time.perf_counter() - start
            # print(f"Elapsed time for resolving sequences: {elapsed:.4f} seconds")
            alignments = self.align_sequences(sequences, references, par=params)
            # elapsed = time.perf_counter() - start
            # print(f"Elapsed time for aligning sequences: {elapsed:.4f} seconds")
            annotated = self.annotate_alignment(
                sequences, references, alignments, par=params
            )
            # elapsed = time.perf_counter() - start
            # print(f"Elapsed time for annotating alignments: {elapsed:.4f} seconds")
            variants = self.annotate_variants(annotated, par=params)
            annotated = self.count_ambigous(annotated, variants)
            write_frames(annotated, alignments_out)
            write_frames(variants, variants_out)
            # raise ValueError()
            return variants

        display = CallSpec(
            path=("BioAligner", "align_and_annotate_variant"),
            kwargs={
                "sequence_paths": sequence_paths,
                "reference_source": self._reference_signature(reference_source),
                "variants_out": variants_out,
                "alignments_out": alignments_out,
                "par": par,
            },
        ).render()
        return Runnable(__call, display=display)

    def count_ambigous(self, annotated: DataFrame, variants: DataFrame, par: ParType | None = None) -> DataFrame:
        annotated = annotated.copy()
        par = self.default_align_params.patch(par)
        query_id_col = par.get("query_id", "R1 Name")
        ambigous_col = par.get("ambigous_col", "is_ambigous")

        number_of_ambigous_variants = variants.groupby(query_id_col)[ambigous_col].sum()
        annotated["n_ambigous_variants"] = annotated[query_id_col].map(number_of_ambigous_variants).fillna(0).astype(int)
        return annotated

    def annotate_ns(self, sequences: DataFrame, par: Params) -> DataFrame:
        query_column = str(par.get("query_sequence", "R1"))
        sequences = sequences.copy()
        sequences["n_count"] = sequences[query_column].str.count("N")
        sequences["contains_n"] = sequences["n_count"] > 0
        return sequences

    def align_sequences(
        self,
        sequences: DataFrame,
        references: DataFrame,
        *,
        par: ParType | None = None,
    ) -> Mapping[str, AlignmentResult]:
        """Align each row in a sequence DataFrame against its reference."""
        par = Params().patch(par)

        sequence_column = str(par.get("query_sequence", "R1"))
        sequence_id_column = str(par.get("query_id", "R1 name"))
        reference_name_column = str(par.get("target_id", "target"))
        reference_sequence_column = str(par.get("target_sequence", "amplicon"))
        if sequence_column not in sequences.columns:
            raise KeyError(
                f"Missing sequence column '{sequence_column}' in input table."
            )
        if reference_name_column not in sequences.columns:
            raise KeyError(
                f"Missing reference name column '{reference_name_column}' in input table."  # noqa: E501
            )
        alignments: Mapping[str, Alignment] = {}
        for idx, row in sequences.iterrows():
            sequence_id = str(row[sequence_id_column])
            query_sequence = self._normalize_sequence(row[sequence_column])
            reference_name = row[reference_name_column]
            reference_sequence = references.loc[
                references[reference_name_column] == reference_name,
                reference_sequence_column,
            ].values[0]
            reference_sequence = self._normalize_sequence(reference_sequence)
            alignment = self._align(reference_sequence, query_sequence)
            alignments[sequence_id] = alignment
        return alignments


    # def annotate_alignments(
    #     self,
    #     sequences: DataFrame,
    #     alignments: Mapping[str, Alignment],
    #     *,
    #     par: ParType | None = None,
    # ) -> DataFrame:

    #     par = self.default_align_params.patch(par)

    #     sequence_column = str(par.get("query_sequence", "sequence"))
    #     sequence_id_column = str(par.get("query_id", "id"))
    #     reference_name_column = str(par.get("target_id", "target"))
    #     chrom_column = str(par.get("chromosome", "chreomosome"))
    #     reference_start_column = str(par.get("start", "start"))
    #     print("annotate")
    #     print(sequence_column, sequence_id_column, reference_name_column, chrom_column, reference_start_column)
    #     if sequence_column not in sequences.columns:
    #         raise KeyError(
    #             f"Missing sequence column '{sequence_column}' in input table."
    #         )
    #     print("sequences.head", sequences.head())
    #     records: list[dict[str, Any]] = []
    #     for idx, row in sequences.iterrows():
    #         reference_name = str(row[reference_name_column])
    #         sequence_id = row[sequence_id_column]
    #         if sequence_id not in alignments:
    #             raise KeyError(f"Missing alignment for sequence_id '{sequence_id}'")
    #         alignment = alignments[sequence_id]

    #         wt_aligned, query_aligned = alignment_to_strings(
    #             alignment,
    #         )

    #         chrom = str(row.get(chrom_column, reference_name))
    #         reference_start = row.get(reference_start_column, 1)

    #         records.append(
    #             {
    #                 "sequence_id": sequence_id,
    #                 "reference_name": reference_name,
    #                 "chrom": chrom,
    #                 "reference_start": reference_start,
    #                 "query": alignment.query,
    #                 "reference": alignment.target,
    #                 "score": float(alignment.score),
    #                 "query_aligned": query_aligned,
    #                 "reference_aligned": wt_aligned,
    #             }
    #         )

    #     return DataFrame.from_records(records)

    def annotate_variants(
        self,
        sequences_aligned: DataFrame,
        *,
        par: ParType | None = None,
    ) -> DataFrame:
        """Extract one-variant-per-row records from an aligned frame."""
        par = self.default_align_params.patch(par)
        variant_rows: list[dict[str, Any]] = []
        chromosome_col = par.get("chromosome", "chromosome")
        reference_start_col = par.get("start", "start")
        query_id_col = par.get("query_id", "R1 Name")
        query_sequence_col = par.get("query_sequence", "R1")
        target_id_col = par.get("target_id", "target")
        target_sequence_col = par.get("target_sequence", "target")
        score_col = par.get("score", "score")
        assembly_col = par.get("assembly", "assembly")
        ambigous_col = par.get("ambigous", "is_ambigous")
        alignment_status_col = "alignment_status"
        for _, row in sequences_aligned.iterrows():
            if row[alignment_status_col] != "success":
                continue
            chrom = row[chromosome_col]
            reference_start = row[reference_start_col]
            assembly = row[assembly_col]
            variants = self.extract_variants_from_gapped(
                reference_aligned=str(row["target_aligned"]),
                query_aligned=str(row["query_aligned"]),
                chrom=chrom,
                reference_start=reference_start,
            )

            for variant in variants:
                variant_row = {
                    query_id_col: row[query_id_col],
                    target_id_col: row[target_id_col],
                    "chrom": variant.chrom,
                    "position": variant.pos,
                    "end": variant.end,
                    "ref": variant.ref,
                    "alt": variant.alt,
                    "type": variant.variant_type,
                    "assembly": assembly,
                    "description": variant.description,
                    score_col: row[score_col],
                    query_sequence_col: row[query_sequence_col],
                    target_sequence_col: row[target_sequence_col],
                    "query_aligned": row["query_aligned"],
                    "reference_aligned": row["target_aligned"],
                    ambigous_col: variant.is_ambiguous,
                }
                variant_rows.append(variant_row)

        if not variant_rows:
            return DataFrame(
                columns=[
                    query_id_col,
                    target_id_col,
                    "chrom",
                    "position",
                    "end",
                    "ref",
                    "alt",
                    "type",
                    "assembly",
                    "description",
                    score_col,
                    query_sequence_col,
                    target_sequence_col,
                    "query_aligned",
                    "reference_aligned",
                ]
            )

        return DataFrame.from_records(variant_rows)

    def extract_variants_from_gapped(
        self,
        *,
        reference_aligned: str,
        query_aligned: str,
        chrom: str,
        reference_start: int = 1,
        reference_sequence: str | None = None,
    ) -> list[Variant]:
        """Extract VCF-conformant variants from aligned reference/query strings."""
        if len(reference_aligned) != len(query_aligned):
            raise ValueError(
                "Aligned strings must have identical length: "
                f"{len(reference_aligned)} != {len(query_aligned)}"
            )
        if reference_sequence is None:
            reference_sequence = reference_aligned.replace("-", "")
        if not reference_sequence:
            raise ValueError("Reference sequence must not be empty.")

        variants: list[Variant] = []
        i = 0
        ref_index = 0

        while i < len(reference_aligned):
            ref_base = reference_aligned[i]
            query_base = query_aligned[i]

            if ref_base == query_base:
                if ref_base == "-":
                    raise ValueError(
                        "Aligned strings must not contain a gap in both rows."
                    )
                ref_index += 1
                i += 1
                continue

            if ref_base == "-" or query_base == "-":
                variant, i, ref_index = self._extract_indel_block(
                    reference_aligned=reference_aligned,
                    query_aligned=query_aligned,
                    chrom=chrom,
                    reference_sequence=reference_sequence,
                    reference_start=reference_start,
                    alignment_index=i,
                    ref_index=ref_index,
                )
                variants.append(variant)
                continue

            variant, i, ref_index = self._extract_substitution_block(
                reference_aligned=reference_aligned,
                query_aligned=query_aligned,
                chrom=chrom,
                reference_sequence=reference_sequence,
                reference_start=reference_start,
                alignment_index=i,
                ref_index=ref_index,
            )
            variants.append(variant)

        return variants

    def _extract_substitution_block(
        self,
        *,
        reference_aligned: str,
        query_aligned: str,
        chrom: str,
        reference_sequence: str,
        reference_start: int,
        alignment_index: int,
        ref_index: int,
    ) -> tuple[Variant, int, int]:
        """Collect directly adjacent substitutions into one SNV or MNV."""
        start = ref_index
        reference_bases: list[str] = []
        query_bases: list[str] = []
        while alignment_index < len(reference_aligned):
            ref_base = reference_aligned[alignment_index]
            query_base = query_aligned[alignment_index]
            if ref_base == "-" or query_base == "-" or ref_base == query_base:
                break
            reference_bases.append(ref_base)
            query_bases.append(query_base)
            alignment_index += 1
            ref_index += 1

        ref = "".join(reference_bases)
        alt = "".join(query_bases)
        pos = reference_start + start
        variant = Variant(
            chrom=chrom,
            pos=pos,
            end=pos + len(ref) - 1,
            ref=ref,
            alt=alt,
            variant_type="SNV" if len(ref) == 1 else "MNV",
            description=f"{pos}{ref}>{alt}",
        )
        self._validate_variant(
            reference_sequence=reference_sequence,
            reference_start=reference_start,
            variant=variant,
        )
        return variant, alignment_index, ref_index

    def _extract_indel_block(
        self,
        *,
        reference_aligned: str,
        query_aligned: str,
        chrom: str,
        reference_sequence: str,
        reference_start: int,
        alignment_index: int,
        ref_index: int,
    ) -> tuple[Variant, int, int]:
        """Collect one contiguous alignment gap and create its VCF indel."""
        if reference_aligned[alignment_index] == "-":
            inserted: list[str] = []
            while alignment_index < len(reference_aligned):
                ref_base = reference_aligned[alignment_index]
                query_base = query_aligned[alignment_index]
                if ref_base != "-" or query_base == "-":
                    break
                inserted.append(query_base)
                alignment_index += 1
            variant = self._make_vcf_indel(
                chrom=chrom,
                reference_sequence=reference_sequence,
                reference_start=reference_start,
                ref_start=ref_index,
                deleted="",
                inserted="".join(inserted),
            )
            return variant, alignment_index, ref_index

        deleted: list[str] = []
        start = ref_index
        while alignment_index < len(reference_aligned):
            ref_base = reference_aligned[alignment_index]
            query_base = query_aligned[alignment_index]
            if ref_base == "-" or query_base != "-":
                break
            deleted.append(ref_base)
            alignment_index += 1
            ref_index += 1
        variant = self._make_vcf_indel(
            chrom=chrom,
            reference_sequence=reference_sequence,
            reference_start=reference_start,
            ref_start=start,
            deleted="".join(deleted),
            inserted="",
        )
        return variant, alignment_index, ref_index

    def _make_vcf_indel(
        self,
        *,
        chrom: str,
        reference_sequence: str,
        reference_start: int,
        ref_start: int,
        deleted: str,
        inserted: str,
    ) -> Variant:
        """Create, left-normalize, and validate an anchored VCF indel."""
        if bool(deleted) == bool(inserted):
            raise ValueError("An indel must have exactly one non-empty allele.")
        if not 0 <= ref_start <= len(reference_sequence):
            raise ValueError(f"Invalid reference index for indel: {ref_start}")

        if ref_start == 0:
            anchor = reference_sequence[0]
            pos = reference_start
            ref = anchor + deleted
            alt = inserted + anchor
        else:
            anchor_index = ref_start - 1
            anchor = reference_sequence[anchor_index]
            pos = reference_start + anchor_index
            ref = anchor + deleted
            alt = anchor + inserted

        pos, ref, alt = self._left_normalize_indel(
            reference_sequence, pos, ref, alt, reference_start
        )
        variant = Variant(
            chrom=chrom,
            pos=pos,
            end=pos + len(ref) - 1,
            ref=ref,
            alt=alt,
            variant_type="DEL" if deleted else "INS",
            description=f"{pos}{ref}>{alt}",
        )
        self._validate_variant(
            reference_sequence=reference_sequence,
            reference_start=reference_start,
            variant=variant,
        )
        return variant

    def _left_normalize_indel(
        self,
        reference_sequence: str,
        pos: int,
        ref: str,
        alt: str,
        reference_start: int = 1,
    ) -> tuple[int, str, str]:
        """Return the minimal leftmost VCF representation of an indel."""
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            pos += 1

        while pos > reference_start and ref[-1] == alt[-1]:
            previous_base = reference_sequence[pos - reference_start - 1]
            ref = previous_base + ref[:-1]
            alt = previous_base + alt[:-1]
            pos -= 1

        return pos, ref, alt

    def _validate_variant(
        self,
        *,
        reference_sequence: str,
        reference_start: int,
        variant: Variant,
    ) -> None:
        """Ensure that the VCF REF allele agrees with the supplied reference."""
        local_pos = variant.pos - reference_start
        expected_ref = reference_sequence[local_pos : local_pos + len(variant.ref)]
        if expected_ref != variant.ref:
            raise ValueError(
                "Variant REF does not match reference sequence: "
                f"expected {expected_ref!r}, got {variant.ref!r} at {variant.pos}"
            )

    def resolve_references(
        self,
        reference_source: Path | ReferenceSpec | Mapping[str, ReferenceSpec],
        par: Params,
    ) -> DataFrame:
        """Resolve references from file, mapping, literal sequence, or spec."""
        par = self.default_align_params.patch(par)
        sequence_col = par.get("target_sequence", "amplicon")
        if isinstance(reference_source, Path):
            frame = read_frame(reference_source)
            self.validate_reference_frame(frame, par=par)
            if sequence_col not in frame.columns:
                api = EnsemblAPI()
                sequences = []
                session = api.open_session()
                for _, row in frame.iterrows():
                    ref_spec = ReferenceSpec(
                        name=row[par.get("target", "target")],
                        species=row[par.get("species", "species")],
                        assembly=row[par.get("assembly", "assembly")],
                        chromosome=row[par.get("chromosome", "chromosome")],
                        start=row[par.get("start", "start")],
                        stop=row[par.get("stop", "stop")],
                        strand=row.get(par.get("strand", "strand"), 1),
                    )
                    sequence = self._normalize_sequence(
                        api.fetch_sequence(ref_spec, session=session)
                    )
                    sequences.append(sequence)
                frame[sequence_col] = sequences
            return frame

        elif isinstance(reference_source, ReferenceSpec):
            api = EnsemblAPI()
            session = api.open_session()
            sequence = api.fetch_sequence(reference_source, session=session)
            frame = DataFrame(
                [
                    {
                        par.get("target", "target"): reference_source.name,  # type: ignore
                        par.get(
                            "species", "species"
                        ): reference_source.species,  # type: ignore
                        par.get(
                            "assembly", "assembly"
                        ): reference_source.assembly,  # type: ignore
                        par.get(
                            "chromosome", "chromosome"
                        ): reference_source.chromosome,  # type: ignore
                        par.get("start", "start"): reference_source.start,  # type: ignore
                        par.get("stop", "stop"): reference_source.stop,  # type: ignore
                        par.get("strand", "strand"): reference_source.strand,  # type: ignore
                        sequence_col: self._normalize_sequence(sequence),
                    }
                ]
            )
            return frame

        elif isinstance(reference_source, Mapping):
            to_frame = []
            for key, value in reference_source.items():
                if isinstance(value, ReferenceSpec):
                    api = EnsemblAPI()
                    session = api.open_session()
                    sequence = api.fetch_sequence(value, session=session)
                    to_frame.append(
                        {
                            par.get("target", "target"): key,
                            par.get("species", "species"): value.species,
                            par.get("assembly", "assembly"): value.assembly,
                            par.get(
                                "chromosome", "chromosome"
                            ): value.chromosome,  # noqa: E501
                            par.get("start", "start"): value.start,
                            par.get("stop", "stop"): value.stop,
                            par.get("strand", "strand"): value.strand,
                            sequence_col: self._normalize_sequence(sequence),
                        }
                    )
            frame = DataFrame(to_frame)
            return frame
        else:
            raise ValueError(
                "Unsupported reference source type. Use Element/FileSource/Path/str/ReferenceSpec/mapping."  # noqa: E501
            )

    def validate_reference_frame(
        self,
        reference_frame: DataFrame,
        par: Params,
    ) -> None:
        for column_name in [
            par.get("target", "target"),
            par.get("species", "species"),
            par.get("assembly", "assembly"),
            par.get("chromosome", "chromosome"),
            par.get("start", "start"),
            par.get("stop", "stop"),
            par.get("strand", "strand"),
        ]:
            if column_name not in reference_frame.columns:
                raise KeyError(f"Missing column '{column_name}' in reference table.")

    def validate_sequence_frame(
        self,
        sequence_frame: DataFrame,
        par: Params,
    ) -> None:
        for column_name in [
            par.get("seq_id", "R1 Name"),
            par.get("target_id", "target"),
            par.get("query", "R1"),
        ]:
            if column_name not in sequence_frame.columns:
                raise KeyError(f"Missing column '{column_name}' in sequence table.")

    def _replace_n_with_reference(
        self,
        query_sequence: str,
        reference_sequence: str,
    ) -> tuple[str, int]:
        """
        Replace ambiguous N bases in the query with the corresponding
        reference base.

        Returns
        -------
        corrected_sequence
            Query sequence with N bases replaced by reference bases.

        n_count
            Number of N bases that were replaced.
        """
        query_sequence = self._normalize_sequence(query_sequence)
        reference_sequence = self._normalize_sequence(reference_sequence)

        if len(query_sequence) != len(reference_sequence):
            raise ValueError(
                "Query and reference must have the same length when "
                "replacing N bases: "
                f"{len(query_sequence)} != {len(reference_sequence)}"
            )

        corrected = []
        n_count = 0

        for query_base, reference_base in zip(
            query_sequence,
            reference_sequence,
        ):
            if query_base.upper() == "N":
                corrected.append(reference_base)
                n_count += 1
            else:
                corrected.append(query_base)

        return "".join(corrected), n_count
    # def _pick_reference_for_row(
    #     self,
    #     *,
    #     row: Any,
    #     references: Mapping[str, str],
    #     reference_name_column: str,
    #     reference_sequence_column: str,
    # ) -> tuple[str, str]:
    #     if (
    #         reference_sequence_column in row
    #         and row[reference_sequence_column] is not None
    #     ):
    #         reference_name = str(row.get(reference_name_column, "inline_reference"))
    #         return (
    #             reference_name,
    #             self._normalize_sequence(row[reference_sequence_column]),
    #         )

    #     row_reference_name = str(row.get(reference_name_column, "default"))
    #     if row_reference_name in references:
    #         return row_reference_name, references[row_reference_name]

    #     if "default" in references:
    #         return row_reference_name, references["default"]

    #     if len(references) == 1:
    #         only_name = next(iter(references.keys()))
    #         return only_name, references[only_name]

    #     raise KeyError(
    #         f"Could not resolve reference for row key '{row_reference_name}'."
    #     )

    # @staticmethod
    # def _source_path(source: Element | FileSource) -> Path:
    #     return source.artifacts.primary.resolve()

    @staticmethod
    def _normalize_sequence(value: Any) -> str:
        seq = str(value).strip().upper()
        if seq.lower() == "nan" or not seq:
            raise ValueError("Encountered empty sequence value.")
        return seq

    # @staticmethod
    # def _to_int(value: Any, *, default: int) -> int:
    #     if value is None:
    #         return default
    #     text = str(value).strip()
    #     if text == "" or text.lower() == "nan":
    #         return default
    #     return int(float(text))

    @staticmethod
    def _reference_signature(
        reference_source: (
            Element
            | FileSource
            | Path
            | str
            | ReferenceSpec
            | Mapping[str, str]
            | Mapping[str, ReferenceSpec]
        ),
    ) -> str:
        if isinstance(reference_source, (Element, FileSource)):
            return str(reference_source.file)
        if isinstance(reference_source, Path):
            return str(reference_source)
        if isinstance(reference_source, str):
            if Path(reference_source).exists():
                return reference_source
        if isinstance(reference_source, ReferenceSpec):
            return (
                f"ReferenceSpec:{reference_source.name}:"  # type: ignore
                f"{len(reference_source.sequence)}"  # type: ignore
            )
            return f"literal_reference:{len(reference_source)}"
        if isinstance(reference_source, Mapping):
            keys = sorted(str(k) for k in reference_source.keys())
            return "mapping:" + ",".join(keys)
        return str(type(reference_source))


def alignment_to_strings(alignment: Alignment) -> tuple[str, str]:
    """Convert a Bio.Align alignment into gapped target/query strings."""

    target, query = alignment.target, alignment.query
    target_blocks, query_blocks = alignment.aligned

    target_aligned = []
    query_aligned = []

    for i, ((t_start, t_end), (q_start, q_end)) in enumerate(
        zip(target_blocks, query_blocks)
    ):
        target_aligned.append(target[t_start:t_end])
        query_aligned.append(query[q_start:q_end])

        if i + 1 == len(target_blocks):
            break

        next_t_start = target_blocks[i + 1][0]
        next_q_start = query_blocks[i + 1][0]

        t_gap = next_t_start - t_end
        q_gap = next_q_start - q_end

        if t_gap > 0 and q_gap == 0:
            target_aligned.append(target[t_end:next_t_start])
            query_aligned.append("-" * t_gap)

        elif t_gap == 0 and q_gap > 0:
            target_aligned.append("-" * q_gap)
            query_aligned.append(query[q_end:next_q_start])

        elif t_gap > 0 and q_gap > 0:
            # Local alignment: skipped region.
            continue

        elif t_gap < 0 or q_gap < 0:
            raise ValueError(
                f"Invalid alignment: target gap={t_gap}, query gap={q_gap}"
            )

    return "".join(target_aligned), "".join(query_aligned)


# @dataclass(frozen=True)
# class Variant:
#     chrom: str
#     pos: int       # 1-based VCF coordinate
#     ref: str
#     alt: str
#     variant_type: str


# class BioAligner:

#     def __init__(self, alignment: str = "pairwise"):
#         self.alignment_type = alignment

#     # ---------------------------------------------------------------------------
# # Alignment
# # ---------------------------------------------------------------------------

#     def create_aligner(
#         self,
#         match_score: float = 2.0,
#         mismatch_score: float = -3.0,
#         gap_open_score: float = -5.0,
#         gap_extend_score: float = -1.0,
#     ) -> PairwiseAligner:
#         """
#         Create a global PairwiseAligner.

#         The parameters are deliberately configurable because the optimal
#         gap penalties depend somewhat on the amplicon/sequence characteristics.
#         """
#         aligner = PairwiseAligner()

#         aligner.mode = "global"

#         aligner.match_score = match_score
#         aligner.mismatch_score = mismatch_score
#         aligner.open_gap_score = gap_open_score
#         aligner.extend_gap_score = gap_extend_score

#         return aligner

#     @element
#     def align2variant(
#         self,
#         sequences: Element | FileSource | FastqSource,
#         reference: Element | FileSource | FastqSource | ReferenceSpec | Mapping[str, ReferenceSpec],
#         *,
#         tag: TagType | None = None,
#         out: OutType | None = None,
#         par: ParType | None = None,
#         cfg: CfgType | None = None,
#     ) -> Element:
#         """
#         Returns an Element that aligns queries to their reference sequence and
#         then extracts the variants from the alignments. the resulting variants
#         are saved in a DataFrame.

#         This columns in the returned frame are:

#         sequence_id query target score query_aligned reference_aligned chrom position alt ref type
#         """
#         # Implement the alignment logic here
#         ref_source = reference
#         if isinstance(reference, (Element, FileSource)):
#             ref_source = reference.file
#         run = self.align_and_annotate_variant(
#             sequences.file,
#             ref_source,
#             par=par,
#         )

#         return Element(
#             ...
#         )


#     def align_and_annotate_variant(
#         self,
#         sequence_path: Path,
#         reference_source: Path | ReferenceSpec | Mapping[str, ReferenceSpec],
#         *,
#         par: ParType | None = None,
#     ) -> Runnable:
#         def run():
#             sequences = read_frame(sequence_path)
#             references = self.resolve_references(par.get("reference", None))()
#             alignments = self.align_sequences(sequences, references, par=par)

#             annotated = self.annotate_alignment(sequences, alignments, par=par)
#             variants = self.annotate_variants(annotated, alignments, par=par)
#             write_frames(variants, out_paths)
#             return annotated

#         callsoec = CallSpec(
#             run=run,
#         )
#         runner = Runnable(run, display=callsoec.render())
#         return runner


#     def align_sequences(
#         self,
#         sequences: DataFrame,
#         references: Mapping[str, str],
#         *,
#         par: ParType | None = None,
#     ) -> Mapping[str, Alignment]:
#         """
#         Aligns queries to their reference sequence and returns a mapping of
#         sequence IDs to their corresponding alignments.

#         Parameters
#         ----------
#         sequences : DataFrame
#             _description_
#         reference : Mapping[str, str]
#             _description_
#         par : ParType | None, optional
#             _description_, by default None

#         Returns
#         -------
#         DataFrame
#             _description_
#         """
#         par = Params().resolve(par)
#         references = self.resolve_references(references)
#         sequence_column = par.get("sequence_column", "sequence")
#         sequence_id_column = par.get("sequence_id_column", "id")
#         scores = []
#         reference_aligned_list = []
#         query_aligned_list = []
#         df = sequences.copy()
#         for _, row in df.iterrows():
#             alignment = self.aligner.align(row["query"], row["reference"])
#             score = alignment.score
#             reference = references[sequence_id]
#             query_aligned, reference_aligned = self.alignment_to_strings(alignment, row["reference"], row["query"])
#             scores.append(score)
#             reference_aligned_list.append(reference_aligned)
#             query_aligned_list.append(query_aligned)

#         df["score"] = scores
#         df["query_aligned"] = query_aligned_list
#         df["reference_aligned"] = reference_aligned_list
#         return df


#    def annotate_alignment(
#         self,
#         sequences: DataFrame,
#         alignments: Mapping[str, Alignment],
#         *,
#         par: ParType | None = None,
#     ) -> DataFrame:
#         """
#         Annotates calculated alignments. The resulting alignments are saved in a
#         DataFrame.

#         This columns in the returned frame are:

#         `sequence_id` `sequence_column` `reference_column` score query_aligned reference_aligned

#         Parameters
#         ----------
#         sequences : DataFrame
#             _description_
#         reference : Mapping[str, str]
#             _description_
#         par : ParType | None, optional
#             _description_, by default None

#         Returns
#         -------
#         DataFrame
#             _description_
#         """
#         par = Params().resolve(par)
#         reference_column = par.get("reference_column", "reference")
#         sequence_column = par.get("sequence_column", "sequence")
#         sequence_id_column = par.get("sequence_id_column", "id")
#         scores = []
#         reference_aligned_list = []
#         query_aligned_list = []
#         df = sequences.copy()
#         for _, row in df.iterrows():
#             sequence_id = row[sequence_id_column]
#             alignment = alignments[sequence_id]
#             score = alignment.score
#             query_aligned, reference_aligned = self.alignment_to_strings(alignment, row["reference"], row["query"])
#             scores.append(score)
#             reference_aligned_list.append(reference_aligned)
#             query_aligned_list.append(query_aligned)

#         df["score"] = scores
#         df["query_aligned"] = query_aligned_list
#         df["reference_aligned"] = reference_aligned_list
#         return df


#     def annotate_variants(
#         self,
#         sequences_aligned: DataFrame,
#         alignments: Mapping[str, Alignment],
#         *,
#         par: ParType | None = None,
#     ) -> Mapping[str, Alignment]:
#         """
#         Annotates calculated alignments. The resulting alignments are saved in a
#         DataFrame.

#         This columns in the returned frame are:

#         `sequence_id` `sequence_column` `reference_column` score query_aligned reference_aligned

#         Parameters
#         ----------
#         sequences : DataFrame
#             _description_
#         reference : Mapping[str, str]
#             _description_
#         par : ParType | None, optional
#             _description_, by default None

#         Returns
#         -------
#         DataFrame
#             _description_
#         """
#         par = Params().resolve(par)
#         reference_column = par.get("reference_column", "reference")
#         sequence_column = par.get("sequence_column", "sequence")
#         sequence_id_column = par.get("sequence_id_column", "id")
#         dfs_to_concat = []
#         for _, row in sequences_aligned.iterrows():
#             sequence_id = row[sequence_id_column]
#             alignment = alignments[sequence_id]
#             variants = self.extract_variants(alignment)
#             for variant in variants:
#                 df_group = row.to_frame().T.copy()
#                 df_group["chrom"] = variant.chrom
#                 df_group["position"] = variant.pos
#                 df_group["ref"] = variant.ref
#                 df_group["alt"] = variant.alt
#                 df_group["type"] = variant.variant_type
#                 dfs_to_concat.append(df_group)
#         return pd.concat(dfs_to_concat, ignore_index=True)


#     def extract_variants(alignment) -> list[Variant]:
#         """
#         Extracts variants from the alignment and returns the Variants as a list.
#         """
#         # Implement variant extraction logic here
#         pass

#     def align_sequence(self, sequence: str, reference: str) -> tuple[str, str]:
#         """
#         Aligns a single sequence to a reference sequence and returns the aligned sequences.
#         """
#         alignment = self.aligner.align(sequence, reference)
#         return alignment, alignment  # Return both query_aligned and reference_aligned as a tuple

#     def resolve_references(
#         self,
#         reference_source: Path | ReferenceSpec | Mapping[str, ReferenceSpec],
#         *,
#         name_column: str | None = None,
#         sequence_column: str | None = None,
#     ) -> Mapping[str, str]:
#         """
#         Resolves the reference sequence from the provided input to a set.
#         """

#         if isinstance(reference_source, ReferenceSpec):
#             # Fetch the reference sequence based on the provided spec
#             sequence = EnsemblAPI.fetch_reference_sequence(reference_source)
#             return defaultdict(lambda: sequence)

#         elif isinstance(reference_source, Path):
#             # Handle other types of references (e.g., Element, FileSource)
#             df = read_frame(reference_source)
#             name_col = name_column if name_column is not None else "name"
#             seq_col = sequence_column if sequence_column is not None else "sequence"
#             return dict(zip(df[name_col], df[seq_col]))

#         elif isinstance(reference_source, Mapping):
#             return reference_source
#         elif isinstance(reference_source, FastqSource):
#             # Handle FastqSource references
#             raise NotImplementedError("FastqSource reference resolution is not implemented yet.")
#         raise ValueError(
#             "Reference must be an Element, FileSource, FastqSource, or ReferenceSpec."
#         )

# def alignment_to_strings(
#     alignment,
#     wt: str,
#     read: str,
# ) -> tuple[str, str]:
#     """
#     Convert a local Bio.Align alignment into gapped WT/read strings.

#     Only gaps that are part of the actual alignment are represented.

#     Regions skipped by the local alignment are NOT interpreted as
#     insertions/deletions.
#     """

#     wt_blocks, read_blocks = alignment.aligned

#     wt_aligned: list[str] = []
#     read_aligned: list[str] = []

#     for i, (
#         (wt_start, wt_end),
#         (read_start, read_end),
#     ) in enumerate(zip(wt_blocks, read_blocks)):

#         # ---------------------------------------------------------
#         # Aligned block
#         # ---------------------------------------------------------

#         wt_block = wt[wt_start:wt_end]
#         read_block = read[read_start:read_end]

#         if len(wt_block) != len(read_block):
#             raise ValueError(
#                 "Aligned block lengths differ: " f"{len(wt_block)} vs {len(read_block)}"
#             )

#         wt_aligned.append(wt_block)
#         read_aligned.append(read_block)

#         # ---------------------------------------------------------
#         # Gap to next alignment block
#         # ---------------------------------------------------------

#         if i + 1 >= len(wt_blocks):
#             continue

#         next_wt_start = wt_blocks[i + 1][0]
#         next_read_start = read_blocks[i + 1][0]

#         wt_gap = next_wt_start - wt_end
#         read_gap = next_read_start - read_end

#         # True deletion:
#         #
#         # WT advances, read does not.
#         #
#         # WT:   ABCDEFG
#         # READ: ABC---G
#         #
#         if wt_gap > 0 and read_gap == 0:

#             wt_aligned.append(wt[wt_end:next_wt_start])

#             read_aligned.append("-" * wt_gap)

#         # True insertion:
#         #
#         # WT:   ABC---DE
#         # READ: ABCXYZDE
#         #
#         elif wt_gap == 0 and read_gap > 0:

#             wt_aligned.append("-" * read_gap)

#             read_aligned.append(read[read_end:next_read_start])

#         # Both sequences advance:
#         #
#         # This is a skipped region caused by the local alignment.
#         #
#         # Do NOT interpret it as an indel.
#         #
#         elif wt_gap > 0 and read_gap > 0:
#             pass

#         elif wt_gap == 0 and read_gap == 0:
#             pass

#         else:
#             raise ValueError(
#                 f"Invalid alignment gap: " f"WT gap={wt_gap}, READ gap={read_gap}"
#             )

#     return (
#         "".join(wt_aligned),
#         "".join(read_aligned),
#     )
