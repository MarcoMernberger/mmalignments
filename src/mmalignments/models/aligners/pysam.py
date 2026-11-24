"""Pysam-backed collection and quantification of blockwise mutations."""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Sequence

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import pysam

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.elements import Element, TableElement, element
from mmalignments.models.overlay import ElementTag
from mmalignments.models.specs import Runnable
from mmalignments.models.tags import Method, Omics, Stage, State
from mmalignments.services.dependencies import depends
from mmalignments.services.io import read_frame, write_frame

_CANDIDATE_COLUMNS = [
    "block",
    "mutation_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "called_samples",
    "filters_by_sample",
]
_COUNT_COLUMNS = [
    "block",
    "mutation_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "sample",
    "called_by_mutect",
    "ref_count",
    "alt_count",
    "other_count",
    "depth",
    "vaf",
    "vaf_percent",
]


def _safe_name(value: str | int) -> str:
    """Return a path-safe identifier while preserving readable block names."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_") or "block"


def _artifact_path(element: Element, key: str) -> Path:
    """Resolve a named element artifact to its concrete file path."""
    artifact = element.artifacts[key]
    if isinstance(artifact, Path):
        return artifact
    path = getattr(artifact, "path", None)
    if isinstance(path, Path):
        return path
    raise TypeError(
        f"Expected the {key!r} artifact of {element.key!r} to be a file path; "
        f"got {type(artifact).__name__}."
    )


def _sample_bam(mutect: Element) -> tuple[str, Path]:
    """Find the single BAM prerequisite belonging to a single-sample call."""
    bams: list[tuple[str, Path]] = []
    for prerequisite in mutect.pres:
        if not isinstance(prerequisite, Element):
            continue
        try:
            bam_path = _artifact_path(prerequisite, "bam")
        except KeyError:
            continue
        bams.append((prerequisite.tag.root, bam_path))

    if len(bams) != 1:
        raise ValueError(
            f"Expected one BAM prerequisite for single-sample Mutect2 element "
            f"{mutect.key!r}, found {len(bams)}."
        )
    return bams[0]


def _validate_samples(
    mutect_elements: Sequence[Element], samples: Sequence[str]
) -> list[Element]:
    """Ensure every requested sample has exactly one matching Mutect2 element."""
    if len(samples) == 0:
        raise ValueError("A block must contain at least one sample.")
    if len(mutect_elements) == 0:
        raise ValueError("At least one Mutect2 element is required.")
    if len(set(samples)) != len(samples):
        raise ValueError("Sample names must be unique within a block.")

    by_sample: dict[str, Element] = {}
    for mutect in mutect_elements:
        sample = mutect.tag.root
        if sample in by_sample:
            raise ValueError(
                f"More than one Mutect2 element was supplied for {sample!r}."
            )
        by_sample[sample] = mutect

    missing = [sample for sample in samples if sample not in by_sample]
    extra = sorted(set(by_sample) - set(samples))
    if missing or extra:
        raise ValueError(
            "Mutect2 elements must match the block's samples exactly; "
            f"missing={missing}, extra={extra}."
        )
    return [by_sample[sample] for sample in samples]


def _read_allele(
    read: pysam.AlignedSegment,
    start: int,
    ref_length: int,
    min_base_quality: int,
) -> tuple[str, float] | None:
    """Reconstruct a read's allele across a VCF reference span from its CIGAR."""
    if read.cigartuples is None or read.query_sequence is None:
        return None

    end = start + ref_length
    reference_position = read.reference_start
    query_position = 0
    covered: set[int] = set()
    allele_parts: list[str] = []
    quality_positions: list[int] = []

    for operation, length in read.cigartuples:
        if operation in (0, 7, 8):  # M, =, X
            overlap_start = max(start, reference_position)
            overlap_end = min(end, reference_position + length)
            if overlap_start < overlap_end:
                query_start = query_position + overlap_start - reference_position
                query_end = query_start + overlap_end - overlap_start
                allele_parts.append(read.query_sequence[query_start:query_end])
                quality_positions.extend(range(query_start, query_end))
                covered.update(range(overlap_start, overlap_end))
            reference_position += length
            query_position += length
        elif operation == 1:  # insertion after the preceding reference base
            if start <= reference_position - 1 < end:
                allele_parts.append(
                    read.query_sequence[query_position : query_position + length]
                )
                quality_positions.extend(range(query_position, query_position + length))
            query_position += length
        elif operation == 2:  # deletion
            overlap_start = max(start, reference_position)
            overlap_end = min(end, reference_position + length)
            if overlap_start < overlap_end:
                covered.update(range(overlap_start, overlap_end))
            reference_position += length
        elif operation == 3:  # reference skip; it is not callable coverage
            reference_position += length
        elif operation == 4:  # soft clip
            query_position += length
        elif operation in (5, 6):  # hard clip or padding
            continue
        else:
            raise ValueError(f"Unsupported CIGAR operation code: {operation}.")

    if len(covered) != ref_length or not quality_positions:
        return None

    qualities = read.query_qualities
    if qualities is None:
        return None
    observed_qualities = [qualities[position] for position in quality_positions]
    if min(observed_qualities) < min_base_quality:
        return None

    mean_quality = sum(observed_qualities) / len(observed_qualities)
    return "".join(allele_parts).upper(), mean_quality


@depends(_read_allele)
def _count_alleles(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    min_mapping_quality: int,
    min_base_quality: int,
    exclude_duplicates: bool,
    ignore_overlapping_mates: bool,
) -> tuple[int, int, int]:
    """Count REF, ALT and other full-span reads for one small variant."""
    if chrom not in bam.references:
        raise ValueError(
            f"Contig {chrom!r} from the Mutect2 VCF is absent from BAM "
            f"{bam.filename!r}."
        )

    observations: list[tuple[str | None, str, float, int, bool]] = []
    start = pos - 1
    for read in bam.fetch(chrom, start, start + len(ref)):
        if (
            read.is_unmapped
            or read.is_secondary
            or read.is_supplementary
            or read.is_qcfail
            or read.mapping_quality == 255
            or read.mapping_quality < min_mapping_quality
            or (exclude_duplicates and read.is_duplicate)
        ):
            continue
        result = _read_allele(read, start, len(ref), min_base_quality)
        if result is None:
            continue
        allele, mean_quality = result
        observations.append(
            (
                read.query_name,
                allele,
                mean_quality,
                read.mapping_quality,
                read.is_read1,
            )
        )

    if ignore_overlapping_mates:
        by_template: dict[str, list[tuple[str | None, str, float, int, bool]]] = (
            defaultdict(list)
        )
        unnamed_observations = []
        for observation in observations:
            query_name = observation[0]
            if query_name is None:
                unnamed_observations.append(observation)
            else:
                by_template[query_name].append(observation)
        observations = unnamed_observations + [
            max(
                template_observations,
                key=lambda item: (
                    item[2],
                    item[3],
                    item[4],
                ),
            )
            for template_observations in by_template.values()
        ]

    ref_count = sum(item[1] == ref.upper() for item in observations)
    alt_count = sum(item[1] == alt.upper() for item in observations)
    other_count = len(observations) - ref_count - alt_count
    return ref_count, alt_count, other_count


class Pysam:
    """Build blockwise Mutect2 candidate and BAM-read-count Elements.

    Candidate variants are the union of records in the supplied, single-sample
    Mutect2 VCFs. Every candidate is then counted in every BAM for the block,
    irrespective of whether that sample called it. Reads failing mapping/base
    quality or alignment-flag filters are excluded. Duplicate-marked reads are
    retained by default, which is appropriate for the current amplicon workflow.
    Overlapping mates are counted once by default, choosing the read with the
    higher mean base quality. Spanning-deletion ``*`` placeholders are ignored;
    other symbolic alleles are rejected because they cannot be represented by
    the small-variant read-counting method.

    VAF is ``ALT / (REF + ALT)``. ``other_count`` is reported separately and
    is not part of that denominator. If no REF or ALT reads pass the filters,
    VAF is missing (NA), not zero. Zero ALT with at least one REF read produces
    a VAF of zero.
    """

    @element
    def collect_mutations(
        self,
        mutect_elements: Sequence[Element],
        *,
        samples: Sequence[str],
        block: str | int,
        outdir: Path = Path("results/mutation_frequencies"),
    ) -> TableElement:
        """Create a per-block union table and record each sample's call status.

        The table contains one row per ``chrom/pos/ref/alt`` allele, the names
        of calling samples, a JSON mapping of each sample's VCF filter labels,
        and one ``called_<sample>`` boolean column per block sample. VCF records
        are intentionally read before Mutect2 filtering so low-confidence
        candidates are not silently lost; their filter labels remain available
        for review.
        """
        ordered_mutects = _validate_samples(mutect_elements, samples)
        sample_names = list(samples)
        columns = _CANDIDATE_COLUMNS + [f"called_{sample}" for sample in sample_names]
        outdir = Path(outdir)
        table_path = outdir / f"{_safe_name(block)}_candidate_variants.tsv"

        def run() -> None:
            calls: dict[tuple[str, int, str, str], dict[str, str]] = defaultdict(dict)
            for sample, mutect in zip(sample_names, ordered_mutects):
                vcf_path = _artifact_path(mutect, "vcf")
                with pysam.VariantFile(str(vcf_path)) as variants:
                    for record in variants:
                        for alt in record.alts or ():
                            if alt == "*":
                                continue
                            if alt.startswith("<"):
                                raise ValueError(
                                    f"Cannot count symbolic allele {alt!r} at "
                                    f"{record.contig}:{record.pos} from {vcf_path}."
                                )
                            key = (
                                record.contig,
                                int(record.pos),
                                record.ref.upper(),
                                alt.upper(),
                            )
                            filters = tuple(record.filter.keys())
                            calls[key][sample] = (
                                ";".join(filters) if filters else "PASS"
                            )

            records: list[dict[str, object]] = []
            for (chrom, pos, ref, alt), filters_by_sample in sorted(calls.items()):
                mutation_id = f"{chrom}:{pos}:{ref}>{alt}"
                row: dict[str, object] = {
                    "block": str(block),
                    "mutation_id": mutation_id,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "called_samples": json.dumps(
                        [
                            sample
                            for sample in sample_names
                            if sample in filters_by_sample
                        ]
                    ),
                    "filters_by_sample": json.dumps(filters_by_sample, sort_keys=True),
                }
                row.update(
                    {
                        f"called_{sample}": sample in filters_by_sample
                        for sample in sample_names
                    }
                )
                records.append(row)

            frame = pd.DataFrame.from_records(records, columns=columns)
            write_frame(frame, table_path, index=False)

        run = depends(_artifact_path, read_frame, write_frame)(run)
        tag = ElementTag(
            root=f"block_{_safe_name(block)}",
            level=1,
            stage=Stage.INTEGRATE,
            method=Method.CUSTOM,
            state=State.MERGED,
            omics=Omics.DNA,
            flag="mutect_candidates",
        )
        inputs = tuple(_artifact_path(mutect, "vcf") for mutect in ordered_mutects)
        return TableElement(
            Element.generate_key(tag, "pysam", "collect_mutations", suffix=str(block)),
            Runnable(run, display=f"pysam collect_mutations(block={block})"),
            tag=tag,
            artifacts=ArtifactSet(table_path, primary_name="tsv"),
            determinants=(f"block={block}", *sample_names),
            inputs=inputs,
            pres=tuple(ordered_mutects),
        )

    @element
    def mutationcount(
        self,
        mutect_elements: Sequence[Element],
        mutation_table: Element,
        *,
        block: str | int,
        min_mapping_quality: int = 20,
        min_base_quality: int = 20,
        exclude_duplicates: bool = False,
        ignore_overlapping_mates: bool = True,
        outdir: Path = Path("results/mutation_frequencies"),
    ) -> TableElement:
        """Count each candidate in every block BAM and write long/wide tables.

        The long ``counts`` table has one row per candidate and sample with
        REF, ALT, other, total depth, VAF, and the original Mutect2 call flag.
        A companion wide VAF matrix has mutations as rows and samples as
        columns. Both outputs preserve missing VAFs for zero evaluable
        REF/ALT reads.
        """
        if min_mapping_quality < 0 or min_base_quality < 0:
            raise ValueError(
                "Mapping and base quality thresholds must be non-negative."
            )
        if not isinstance(mutation_table, Element):
            raise TypeError("mutation_table must be an Element.")

        sample_names = [mutect.tag.root for mutect in mutect_elements]
        ordered_mutects = _validate_samples(mutect_elements, sample_names)
        table_path = _artifact_path(mutation_table, "tsv")
        required = set(
            _CANDIDATE_COLUMNS + [f"called_{sample}" for sample in sample_names]
        )
        outdir = Path(outdir)
        stem = f"block_{_safe_name(block)}"
        counts_path = outdir / f"{stem}_mutation_counts.tsv"
        matrix_path = outdir / f"{stem}_vaf_matrix.tsv"
        bam_paths = tuple(_sample_bam(mutect)[1] for mutect in ordered_mutects)

        def run() -> None:
            candidates = read_frame(table_path)
            missing_columns = sorted(required - set(candidates.columns))
            if missing_columns:
                raise ValueError(
                    f"Candidate table {table_path} is missing columns: "
                    f"{missing_columns}."
                )
            if not candidates.empty and set(candidates["block"].astype(str)) != {
                str(block)
            }:
                raise ValueError(
                    f"Candidate table {table_path} does not contain only block "
                    f"{block!r} variants."
                )
            rows: list[dict[str, object]] = []
            for sample, bam_path in zip(sample_names, bam_paths):
                with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                    for candidate in candidates.to_dict(orient="records"):
                        ref_count, alt_count, other_count = _count_alleles(
                            bam,
                            str(candidate["chrom"]),
                            int(candidate["pos"]),
                            str(candidate["ref"]),
                            str(candidate["alt"]),
                            min_mapping_quality=min_mapping_quality,
                            min_base_quality=min_base_quality,
                            exclude_duplicates=exclude_duplicates,
                            ignore_overlapping_mates=ignore_overlapping_mates,
                        )
                        evaluable_depth = ref_count + alt_count
                        vaf = (
                            alt_count / evaluable_depth
                            if evaluable_depth
                            else float("nan")
                        )
                        rows.append(
                            {
                                "block": str(block),
                                "mutation_id": candidate["mutation_id"],
                                "chrom": candidate["chrom"],
                                "pos": int(candidate["pos"]),
                                "ref": candidate["ref"],
                                "alt": candidate["alt"],
                                "sample": sample,
                                "called_by_mutect": candidate[f"called_{sample}"],
                                "ref_count": ref_count,
                                "alt_count": alt_count,
                                "other_count": other_count,
                                "depth": ref_count + alt_count + other_count,
                                "vaf": vaf,
                                "vaf_percent": vaf * 100,
                            }
                        )

            counts = pd.DataFrame.from_records(rows, columns=_COUNT_COLUMNS)
            if not counts.empty:
                counts["vaf_percent"] = counts["vaf"] * 100
            vaf_matrix = counts.pivot(
                index="mutation_id", columns="sample", values="vaf_percent"
            )
            vaf_matrix.index.name = "mutation_id"
            vaf_matrix = vaf_matrix.reindex(columns=sample_names)
            write_frame(counts, counts_path, index=False)
            write_frame(vaf_matrix, matrix_path, index=True)

        run = depends(_count_alleles, _sample_bam, read_frame, write_frame)(run)
        tag = ElementTag(
            root=f"block_{_safe_name(block)}",
            level=2,
            stage=Stage.QUANT,
            method=Method.CUSTOM,
            state=State.COUNT,
            omics=Omics.DNA,
            flag="mutation_frequencies",
        )
        return TableElement(
            Element.generate_key(
                tag,
                "pysam",
                "mutationcount",
                suffix=(
                    f"{block}_mapq{min_mapping_quality}_baseq{min_base_quality}_"
                    f"dedup{exclude_duplicates}_overlap{ignore_overlapping_mates}"
                ),
            ),
            Runnable(run, display=f"pysam mutationcount(block={block})"),
            tag=tag,
            artifacts=ArtifactSet(counts_path, primary_name="tsv").with_extra(
                "vaf_matrix", matrix_path
            ),
            determinants=(
                f"min_mapping_quality={min_mapping_quality}",
                f"min_base_quality={min_base_quality}",
                f"exclude_duplicates={exclude_duplicates}",
                f"ignore_overlapping_mates={ignore_overlapping_mates}",
            ),
            inputs=(table_path, *bam_paths),
            pres=(mutation_table, *ordered_mutects),
        )

    @element
    def plot_mutation_frequencies(
        self,
        counts: Element,
        *,
        block: str | int,
        outdir: Path = Path("results/mutation_frequencies"),
    ) -> Element:
        """Create an interactive grouped bar plot of per-sample VAF percentages.

        Bars are grouped by mutation and colored by sample. Missing VAF values
        (no evaluable REF/ALT reads) are omitted by Plotly; a measured zero VAF
        remains represented as a zero-height bar.
        """
        counts_path = counts.artifacts["counts"].resolve()
        required = {"mutation_id", "sample", "vaf_percent"}
        outdir = Path(outdir)
        plot_path = outdir / f"{_safe_name(block)}_vaf_barplot.html"

        def run() -> None:
            frame = read_frame(counts_path)
            missing_columns = sorted(required - set(frame.columns))
            if missing_columns:
                raise ValueError(
                    f"Counts table {counts_path} is missing columns: "
                    f"{missing_columns}."
                )
            if frame.empty:
                figure = go.Figure()
                figure.add_annotation(
                    text="No Mutect2 candidate variants were found for this block.",
                    showarrow=False,
                )
            else:
                figure = px.bar(
                    frame,
                    x="mutation_id",
                    y="vaf_percent",
                    color="sample",
                    barmode="group",
                    hover_data=[
                        column
                        for column in (
                            "chrom",
                            "pos",
                            "ref",
                            "alt",
                            "ref_count",
                            "alt_count",
                            "other_count",
                            "depth",
                            "called_by_mutect",
                        )
                        if column in frame.columns
                    ],
                    labels={
                        "mutation_id": "Mutation",
                        "vaf_percent": "VAF (%)",
                        "sample": "Sample",
                    },
                    title=f"Block {block}: mutation frequencies across samples",
                )
            figure.update_layout(
                xaxis_title="Mutation (chromosome:position:REF>ALT)",
                yaxis_title="VAF (%)",
                yaxis_range=[0, 100],
                legend_title="Sample",
            )
            plot_path.parent.mkdir(parents=True, exist_ok=True)
            figure.write_html(plot_path, include_plotlyjs=True)

        run = depends(read_frame)(run)
        tag = ElementTag(
            root=f"block_{_safe_name(block)}",
            level=3,
            stage=Stage.SUMMARY,
            method=Method.CUSTOM,
            state=State.PLOT,
            omics=Omics.DNA,
            flag="mutation_frequencies",
        )
        return Element(
            Element.generate_key(
                tag, "pysam", "plot_mutation_frequencies", suffix=str(block)
            ),
            Runnable(
                run,
                display=f"pysam plot_mutation_frequencies(block={block})",
            ),
            tag=tag,
            artifacts=ArtifactSet(plot_path, primary_name="html"),
            inputs=(counts_path,),
            pres=(counts,),
        )
