"""Wrapper for Ensembl VEP in the mmalignments Element design.

This module provides two main entry points:

- :meth:`Vep.build_input`: create a VEP-compatible input file from a variant
  table (or from mutated sequences + optional WT sequence)
- :meth:`Vep.run`: run the ``vep`` CLI for one prepared input file
"""

from __future__ import annotations

import logging
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, cast

import pandas as pd  # type: ignore[import]
import pysam  # type: ignore[import]
from pandas import DataFrame, Series  # type: ignore[import]

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.elements import CallSpec, Element, FileSource, element
from mmalignments.models.overlay import (
    CfgType,
    FileSpec,
    Out,
    OutSpec,
    OutType,
    Params,
    ParType,
    TagType,
)
from mmalignments.models.parameters import ParamSet
from mmalignments.models.tags import Method, Stage, State
from mmalignments.services.io import parents, read_frame, write_frame

from ..externals import External, Runnable, SubroutineIn, subroutine

logger = logging.getLogger(__name__)


def _params_to_cli(
    params: ParType | None,
    *,
    skip_keys: set[str] | None = None,
) -> list[str]:
    """Convert overlay params to CLI flags.

    Keys are converted from ``snake_case`` to ``--kebab-case``.
    ``True`` values are emitted as flag-only options.
    ``False``/``None``/empty lists are skipped.
    """
    if params is None:
        return []

    skip = skip_keys or set()
    cli: list[str] = []
    for key in sorted(params.keys()):
        if key in skip:
            continue
        value = params.get(key)
        if value is None or value is False or value == []:
            continue
        flag = "--" + str(key).replace("_", "-")
        if value is True:
            cli.append(flag)
            continue
        if isinstance(value, (list, tuple)):
            cli.extend([flag, ",".join(str(v) for v in value)])
            continue
        cli.extend([flag, str(value)])
    return cli


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    text = str(value).strip().lower()
    return text in {"", "none", "null"}


def _clean_string(value: Any, default: str = "") -> str:
    return default if _is_missing(value) else str(value).strip()


def _parse_haplotype_variant_row(
    row: Any,
    *,
    sample_column: str,
    chrom_column: str,
    position_column: str,
    ref_column: str,
    alt_column: str,
    ambigous_column: str = "is_ambiguous",
) -> tuple[str, tuple[str, int, str, str], bool]:
    """Validate and turn one variant-table row into a haplotype allele key."""
    sample = _clean_string(row[sample_column])
    chrom = _clean_string(row[chrom_column])
    ref = _clean_string(row[ref_column]).upper()
    alt = _clean_string(row[alt_column]).upper()
    is_ambiguous = bool(row[ambigous_column])
    if not sample or not chrom or not ref or not alt:
        raise ValueError(
            f"Variant rows require sample, chrom, ref, and alt, was {row}."
        )
    if "-" in ref or "-" in alt:
        raise ValueError("VCF alleles must not contain alignment gaps.")
    position = int(row[position_column])
    if position < 1:
        raise ValueError(f"VCF position must be positive: {position}")
    return sample, (chrom, position, ref, alt), is_ambiguous


def _write_haplotype_vcf(
    output_vcf: Path,
    sample_variants: dict[str, set[tuple[str, int, str, str]]],
    assembly: str,
) -> None:
    """Write a bgzip-compressed multi-sample VCF of haploid haplotypes."""
    sorted_variants = sorted(
        set().union(*sample_variants.values()) if sample_variants else set()
    )
    header = pysam.VariantHeader()
    # header.add_meta("fileformat", value="VCFv4.3")
    header.add_meta("source", value="mmalignments")
    header.add_meta("reference", value=assembly)
    header.add_meta("phasing", value="partial")
    for contig in sorted({variant[0] for variant in sorted_variants}):
        header.contigs.add(contig)
    header.add_meta(
        "INFO",
        items=[
            ("ID", "TYPE"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Variant type: SNV/MNV/INS/DEL"),
        ],
    )
    header.add_meta(
        "FORMAT",
        items=[
            ("ID", "GT"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Haploid genotype"),
        ],
    )
    for sample in sorted(sample_variants):
        header.add_sample(sample)

    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    with pysam.VariantFile(str(output_vcf), "wz", header=header) as vcf:
        for chrom, position, ref, alt in sorted_variants:
            record = vcf.new_record(
                contig=chrom,
                start=position - 1,
                stop=position - 1 + len(ref),
                alleles=(ref, alt),
            )
            record.info["TYPE"] = variant_type_from_alleles(ref, alt)
            variant_key = (chrom, position, ref, alt)
            for sample, variants in sample_variants.items():
                genotype = 1 if variant_key in variants else 0
                record.samples[sample]["GT"] = (genotype, genotype)
                record.samples[sample].phased = True
            vcf.write(record)
    pysam.tabix_index(str(output_vcf), preset="vcf", force=True)


def _to_strand_token(strand: str) -> str:
    return "+" if str(strand).strip() not in {"-", "-1"} else "-"


def _to_region_strand(strand: str) -> str:
    return "-1" if _to_strand_token(strand) == "-" else "1"


@dataclass(frozen=True)
class _SimpleVariant:
    chrom: str
    start: int
    end: int
    ref: str
    alt: str
    strand: str
    variant_id: str | None = None


def _diff_variant_from_sequences(
    wt_sequence: str,
    alt_sequence: str,
    genomic_start: int,
) -> _SimpleVariant | None:
    """Derive a single minimal variant from WT and altered sequence.

    ``genomic_start`` is expected to be 1-based coordinate of the first base
    in the sequence window.
    """
    wt = wt_sequence.upper()
    alt = alt_sequence.upper()
    if wt == alt:
        return None

    prefix = 0
    min_len = min(len(wt), len(alt))
    while prefix < min_len and wt[prefix] == alt[prefix]:
        prefix += 1

    suffix = 0
    wt_remaining = len(wt) - prefix
    alt_remaining = len(alt) - prefix
    while (
        suffix < wt_remaining
        and suffix < alt_remaining
        and wt[len(wt) - 1 - suffix] == alt[len(alt) - 1 - suffix]
    ):
        suffix += 1

    wt_mid = wt[prefix : len(wt) - suffix]
    alt_mid = alt[prefix : len(alt) - suffix]
    pos = genomic_start + prefix

    if wt_mid == "":
        # Ensembl format insertion convention: start = end + 1.
        return _SimpleVariant(
            chrom="",
            start=pos,
            end=pos - 1,
            ref="-",
            alt=alt_mid,
            strand="+",
        )
    if alt_mid == "":
        return _SimpleVariant(
            chrom="",
            start=pos,
            end=pos + len(wt_mid) - 1,
            ref=wt_mid,
            alt="-",
            strand="+",
        )

    return _SimpleVariant(
        chrom="",
        start=pos,
        end=pos + len(wt_mid) - 1,
        ref=wt_mid,
        alt=alt_mid,
        strand="+",
    )


def _to_vcf_ref_alt(
    variant: _SimpleVariant,
    wt_sequence: str | None = None,
    genomic_start: int | None = None,
) -> tuple[int, str, str]:
    """Convert a simple variant to POS/REF/ALT for VCF.

    VCF needs anchored alleles for insertions/deletions. If sequence context is
    missing, a conservative ``N`` anchor is used.
    """
    if variant.ref != "-" and variant.alt != "-":
        return variant.start, variant.ref, variant.alt

    anchor_base = "N"
    if wt_sequence is not None and genomic_start is not None:
        idx = variant.start - genomic_start
        if variant.ref == "-":
            # insertion between idx-1 and idx
            if 0 < idx <= len(wt_sequence):
                anchor_base = wt_sequence[idx - 1]
            elif 0 <= idx < len(wt_sequence):
                anchor_base = wt_sequence[idx]
        else:
            # deletion starts at idx
            if idx > 0 and (idx - 1) < len(wt_sequence):
                anchor_base = wt_sequence[idx - 1]
            elif 0 <= idx < len(wt_sequence):
                anchor_base = wt_sequence[idx]

    if variant.ref == "-":
        pos = max(1, variant.start - 1)
        ref = anchor_base
        alt = anchor_base + variant.alt
        return pos, ref, alt

    pos = max(1, variant.start - 1)
    ref = anchor_base + variant.ref
    alt = anchor_base
    return pos, ref, alt


def _render_line(
    variant: _SimpleVariant,
    output_format: str,
    *,
    wt_sequence: str | None = None,
    genomic_start: int | None = None,
) -> str:
    """Render a single variant line for one VEP input format."""
    fmt = output_format.lower()

    if fmt in {"ensembl", "default"}:
        parts = [
            variant.chrom,
            str(variant.start),
            str(variant.end),
            f"{variant.ref}/{variant.alt}",
            _to_strand_token(variant.strand),
        ]
        if variant.variant_id:
            parts.append(variant.variant_id)
        return "\t".join(parts)

    if fmt == "vcf":
        pos, ref, alt = _to_vcf_ref_alt(
            variant,
            wt_sequence=wt_sequence,
            genomic_start=genomic_start,
        )
        vcf_id = variant.variant_id or "."
        return "\t".join([variant.chrom, str(pos), vcf_id, ref, alt, ".", ".", "."])

    if fmt in {"id", "identifier"}:
        if not variant.variant_id:
            raise ValueError("id format requires a variant identifier per row")
        return variant.variant_id

    if fmt == "region":
        allele = variant.alt if variant.alt != "-" else "DEL"
        return (
            f"{variant.chrom}:{variant.start}-{variant.end}:"
            f"{_to_region_strand(variant.strand)}/{allele}"
        )

    raise ValueError(f"Unsupported VEP input format: {output_format}")


class Vep(External):
    """Wrapper for the ``vep`` CLI."""

    def __init__(
        self,
        name: str = "vep",
        primary_binary: str = "vep",
        version: str | None = None,
        source: str = "https://www.ensembl.org/info/docs/tools/vep/",
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
    ) -> None:
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=parameters or {},
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        """Detect VEP version from ``vep --help`` output."""
        if self._version:
            return self._version
        if not self.primary_binary or not self.ensure_binary():
            return fallback
        try:
            cp = subprocess.run(
                [self.primary_binary, "--help"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            out = (cp.stdout or "") + "\n" + (cp.stderr or "")
            m = re.search(r"ensembl-vep\s*:\s*([0-9]+(?:\.[0-9]+)*)", out)
            if m:
                return m.group(1)
        except Exception:
            pass
        return fallback

    def default_outdir(self, sample_name: str) -> Path:
        """Default result location for VEP output files."""
        return Path("results") / "vep" / self.version_name / sample_name

    @property
    def vep_params(self) -> Params:
        return Params()

    @property
    def haplo_params(self) -> Params:
        return Params()

    @element
    def variants2vcf(
        self,
        variants: Element | FileSource,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Create a vcf file from variants dataframe"""
        tag = variants.tag.bump(
            stage=Stage.ANALYSIS,
            method=Method.VEP,
            state=State.PROCESSED,
            flag="variant2vcf",
        ).resolve(tag)
        out = OutSpec.from_tag(
            tag,
            Out(
                folder=self.default_outdir(variants.tag.root),
                ext="vcf.gz",
            ),
        )
        out.add_output("tbx", FileSpec(out.file.name, ext="vcf.gz.tbi"))
        artifacts = ArtifactSet.from_outspec(out)
        params = Params().resolve(par)  # type: ignore
        runner = self.write_variants_to_vcf(
            source_path=variants.file,
            out_paths=out.files,
            par=params,
        )

        key = Element.generate_key(tag, "VEP", "variant2vcf")
        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=params.determinants,
            pres=(variants,),
        )

    def write_variants_to_vcf(
        self,
        source_path: Path,
        out_paths: list[Path],
        omit_ambigous: bool = True,
        *,
        par: ParType | None = None,
    ) -> Runnable:
        """Write one haploid, phased VCF sample per original query sequence."""

        def run():
            frame = read_frame(source_path)
            params = Params.resolve(par)  # type: ignore
            sample_column = str(params.get("sample_column", "R1 Name"))
            chrom_column = str(params.get("chrom_column", "chrom"))
            position_column = str(params.get("position_column", "position"))
            ref_column = str(params.get("ref_column", "ref"))
            alt_column = str(params.get("alt_column", "alt"))
            assembly_column = str(params.get("assembly_column", "assembly"))
            ambigous_column = str(params.get("ambigous_column", "is_ambigous"))

            required_columns = {
                sample_column,
                chrom_column,
                position_column,
                ref_column,
                alt_column,
                assembly_column,
                ambigous_column,
            }
            missing_columns = required_columns.difference(frame.columns)
            if missing_columns:
                raise KeyError(
                    "Variant table is missing required columns: "
                    + ", ".join(sorted(missing_columns))
                )
            if not out_paths:
                raise ValueError("At least one output path is required for the VCF.")
            assemblies = {
                _clean_string(assembly)
                for assembly in frame[assembly_column]
                if not _is_missing(assembly)
            }
            if len(assemblies) != 1:
                raise ValueError(
                    "Variant table must contain exactly one non-empty assembly."
                )
            assembly = assemblies.pop()

            sample_variants: dict[str, set[tuple[str, int, str, str]]] = {}
            for _, row in frame.iterrows():
                sample, variant_key, is_ambiguous = _parse_haplotype_variant_row(
                    row,
                    sample_column=sample_column,
                    chrom_column=chrom_column,
                    position_column=position_column,
                    ref_column=ref_column,
                    alt_column=alt_column,
                    ambigous_column=ambigous_column,
                )
                if omit_ambigous and is_ambiguous:
                    continue
                sample_variants.setdefault(sample, set()).add(variant_key)
            _write_haplotype_vcf(out_paths[0], sample_variants, assembly)

        callspec = CallSpec(
            path=("write_variants_to_vcf",),
            kwargs={"source_path": source_path, "out_paths": out_paths, "par": par},
        )
        return Runnable(
            run,
            display=callspec.render(),
        )

    @element
    def vep(
        self,
        vep_input: Element | FileSource,
        *,
        input_format: str = "vcf",
        species: str = "mus_musculus",
        transcript_id: str | None = None,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Run VEP for one prepared input file."""
        par = self.vep_params.resolve(par)  # type: ignore
        tag = vep_input.tag.bump(
            stage=Stage.ANALYSIS, method=Method.VEP, state=State.ANNOTATED, flag=species
        ).resolve(tag)
        out = OutSpec.from_tag(
            tag,
            folder=self.default_outdir(tag.root),
            ext="txt",
        ).resolve(out)
        artifacts = ArtifactSet.from_outspec(out)
        runner = self.run_vep(
            input_file=vep_input.file,
            output_file=out.file,
            input_format=input_format,
            species=species,
            transcript_id=transcript_id,
            par=par,
            cfg=cfg,
        )
        key = Element.generate_key(
            tag,
            "run",
            "vep",
        )
        determinants = (species, input_format) + par.determinants

        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=(vep_input.file,),
            pres=(vep_input,),
        )

    @subroutine
    def run_vep(
        self,
        input_file: Path,
        output_file: Path,
        *,
        input_format: str = "vcf",
        species: str = "homo_sapiens",
        transcript_id: str | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for one ``vep`` invocation."""
        _ = cfg  # consumed by @subroutine wrapper
        par = self.vep_params.resolve(par)  # type: ignore

        vep_format = "ensembl" if input_format == "default" else input_format
        arguments = [
            "--input_file",
            str(Path(input_file).resolve()),
            "--output_file",
            str(Path(output_file).resolve()),
            "--species",
            str(species),
            "--format",
            str(vep_format),
            "--force_overwrite",
            "--database",
            "--hgvs",
            "--hgvsg",
        ]
        if transcript_id is not None:
            arguments.extend(["--transcript_filter", str(transcript_id)])
        arguments.extend(
            _params_to_cli(
                par,
            )
        )

        in_paths = [Path(input_file).resolve()]
        out_paths = [Path(output_file).resolve()]
        return (
            arguments,
            None,
            in_paths,
            out_paths,
            None,
            None,
            None,
        )

    @element
    def annotate(
        self,
        variant: Element,
        haplo: Element,
        vep: Element,
        *,
        transcript: str | None = None,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Run VEP for one prepared input file."""
        par = self.vep_params.resolve(par)  # type: ignore
        tag = variant.tag.bump(
            stage=Stage.ANALYSIS,
            method=Method.VEP,
            state=State.ANNOTATED,
            flag="vep.haplo",
        ).resolve(tag)
        out = OutSpec.from_tag(
            tag,
            folder=self.default_outdir(tag.root),
            ext="tsv",
        ).resolve(out)
        artifacts = ArtifactSet.from_outspec(out)
        runner = self.combine_vep_haplo(
            sequence_path=variant.artifacts["alignments"].resolve(),
            vep_path=vep.file,
            haplo_path=haplo.file,
            variant_path=variant.file,
            output_path=out.file,
            transcript=transcript,
            par=par,
        )
        key = Element.generate_key(
            tag,
            "vep",
            "annotate",
        )
        determinants = par.determinants
        if transcript is not None:
            determinants += (transcript,)

        return Element(
            key,
            run=runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=(),
            pres=(variant, vep, haplo),
        )

    def combine_vep_haplo(
        self,
        sequence_path: Path,
        vep_path: Path,
        haplo_path: Path,
        variant_path: Path,
        output_path: Path,
        *,
        transcript: str | None = None,
        rename_dict: dict[str, str] | None = None,
        par: ParType | None = None,
    ) -> Runnable:
        """
        Combine sequence, VEP and Haplosaurus annotations.

        Haplosaurus annotations are always added when available, regardless
        of the number of variants in the original sequence.

        Parameters
        ----------
        sequence_path
            Sequence dataframe containing exactly one row per original
            sequence.

        vep_path
            VEP annotation output.

        haplo_path
            Haplosaurus annotation output.

        variant_path
            Variant dataframe containing one row per variant.

        output_path
            Output dataframe.

        transcript
            Optional Ensembl transcript ID used to filter VEP and Haplosaurus
            annotations.

        rename_dict
            Optional mapping used to replace Ensembl IDs with manually supplied
            identifiers, e.g.

                {
                    "ENSMUST00000108658": "NM_011640.3",
                    "ENSMUSP00000104298": "NP_035770.2",
                }

            The dictionary is applied to the final VEP/Haplosaurus annotation
            strings only. VEP/Haplosaurus lookup continues to use Ensembl IDs.

        par
            Optional parameter override.

        Returns
        -------
        Runnable
            Runnable producing the combined annotation dataframe.

        Notes
        -----
        ``vep_status`` and ``haplo_status`` describe whether an annotation was
        successfully found.

        Possible statuses:

            VEP:
                "annotated"
                "not_found"
                "empty"

            Haplosaurus:
                "annotated"
                "not_found"
                "empty"

        ``annotation_source`` is one of:

            "VEP"
            "Haplosaurus"
            "VEP+Haplosaurus"
            "none"
        """

        param = self.vep_params.resolve(par)  # type: ignore
        sample_column = cast(str, param.get("sample_column", "R1 Name"))
        chrom_column = cast(str, param.get("chrom_column", "chrom"))
        position_column = cast(str, param.get("position_column", "position"))
        ref_column = cast(str, param.get("ref_column", "ref"))
        alt_column = cast(str, param.get("alt_column", "alt"))
        ambiguous_column = cast(
            str, param.get("ambiguous_column", "n_ambigous_variants")
        )

        if transcript is not None:
            transcript = str(transcript)

        def run(
            transcript=transcript,
            rename_dict=rename_dict,
        ):
            # Read input
            sequence_frame = read_frame(sequence_path)
            variant_frame = read_frame(variant_path)
            vep_frame = _read_vep_output(vep_path)
            haplo_frame = _read_haplo_output(haplo_path)
            # Validate sequence dataframe
            if sample_column not in sequence_frame.columns:
                raise KeyError(
                    f"Sequence dataframe is missing sample column "
                    f"{sample_column!r}."
                )

            # Validate variant dataframe
            required_variant_columns = {
                sample_column,
                chrom_column,
                position_column,
                ref_column,
                alt_column,
            }

            missing = required_variant_columns.difference(variant_frame.columns)

            if missing:
                raise KeyError(
                    "Variant dataframe is missing required columns: "
                    + ", ".join(sorted(missing))
                )

            # Group variants by original sequence
            sequence_variants = build_sequence_variants(
                variant_frame,
                sample_column,
                chrom_column,
                position_column,
                ref_column,
                alt_column,
            )

            # Build VEP/haplo lookup
            vep_lookup = build_vep_lookup(vep_frame, transcript=transcript)
            haplo_lookup = build_haplosaurus_lookup(haplo_frame, transcript=transcript)

            # Create one annotation row per sequence
            annotation_rows = build_annotation_rows(
                sequence_frame,
                sequence_variants,
                sample_column,
                vep_lookup,
                haplo_lookup,
                rename_dict=rename_dict,
                ambiguous_column=ambiguous_column,
            )
            annotation_frame = DataFrame(annotation_rows)

            # Merge annotation onto sequence dataframe.
            result = sequence_frame.merge(
                annotation_frame,
                on=sample_column,
                how="left",
                validate="one_to_one",
            )

            # Write result
            parents(output_path)
            write_frame(
                result,
                output_path,
            )

        callspec = CallSpec(
            path=("combine_vep_haplo",),
            kwargs={
                "sequence_path": sequence_path,
                "variant_path": variant_path,
                "vep_path": vep_path,
                "haplo_path": haplo_path,
                "output_path": output_path,
                "par": par,
            },
        )

        return Runnable(
            run,
            display=callspec.render(),
        )


def build_sequence_variants(
    variant_frame: DataFrame,
    sample_column: str,
    chrom_column: str,
    position_column: str,
    ref_column: str,
    alt_column: str,
) -> dict[str, list[tuple[str, int, str, str]]]:
    """
    build_sequence_variants builds a dictionary of sequence variants keyed by sample ID.

    Parameters
    ----------
    variant_frame : DataFrame
        DataFrame containing variant information.
    sample_column : str
        Name of the column containing sample IDs.
    chrom_column : str
        Name of the column containing chromosome information.
    position_column : str
        Name of the column containing variant positions.
    ref_column : str
        Name of the column containing reference alleles.
    alt_column : str
        Name of the column containing alternate alleles.

    Returns
    -------
    dict[str, list[tuple[str, int, str, str]]]
        Dictionary keyed by sample ID, with values being lists of tuples representing
        variants. Each tuple contains (chromosome, position, reference allele,
        alternate allele).

    Raises
    ------
    ValueError
        Raised if a variant row contains an empty sample ID.
    """
    sequence_variants: dict[
        str,
        list[tuple[str, int, str, str]],
    ] = {}

    for _, row in variant_frame.iterrows():
        sample = _clean_string(row[sample_column])

        if not sample:
            raise ValueError("Variant row contains an empty sample ID.")

        variant = _variant_key_from_row(
            row,
            chrom_column=chrom_column,
            position_column=position_column,
            ref_column=ref_column,
            alt_column=alt_column,
        )

        sequence_variants.setdefault(
            sample,
            [],
        ).append(variant)

    # Remove duplicate variants while preserving order.
    for sample, variants in sequence_variants.items():
        sequence_variants[sample] = list(dict.fromkeys(variants))
    return sequence_variants


def build_vep_lookup(
    vep_frame: DataFrame, transcript: str | None = None
) -> dict[str, dict[str, str]]:
    """
    build_vep_lookup builds a lookup dictionary for VEP annotations keyed by uploaded
    variation.

    This function processes a VEP DataFrame and constructs a nested dictionary where the
    first level keys are the uploaded variation identifiers and the second level keys
    are annotation types (hgvs_g, hgvs_c, hgvs_p) with their corresponding values.

    Parameters
    ----------
    vep_frame : DataFrame
        DataFrame containing VEP annotations.
    transcript : str | None, optional
        Transcript ID to filter the VEP annotations, by default None.

    Returns
    -------
    dict[str, dict[str, str]]
        Dictionary keyed by uploaded variation, with values being dictionaries of VEP
        annotations.

    Raises
    ------
    KeyError
        Raised if the VEP DataFrame is missing required columns.
    """
    required_vep_columns = {
        "Uploaded_variation",
        "Feature",
        "Extra",
    }

    missing = required_vep_columns.difference(vep_frame.columns)

    if missing:
        raise KeyError(
            "VEP output is missing required columns: " + ", ".join(sorted(missing))
        )

    if transcript is not None:
        vep_frame = vep_frame.loc[
            vep_frame["Feature"].eq(transcript)
        ].copy()  # type: ignore

    vep_lookup: dict[
        str,
        dict[str, str],
    ] = {}

    for _, row in vep_frame.iterrows():
        uploaded_variation = _clean_string(row["Uploaded_variation"])

        extra = _parse_vep_extra(row["Extra"])

        vep_annotation = {
            "hgvs_g": _clean_string(extra.get("HGVSg", "")),
            "hgvs_c": _clean_string(extra.get("HGVSc", "")),
            "hgvs_p": _clean_string(extra.get("HGVSp", "")),
        }

        # If the same uploaded variation occurs more than once,
        # preserve the first annotation as before.
        vep_lookup.setdefault(
            uploaded_variation,
            vep_annotation,
        )
    return vep_lookup


def build_haplosaurus_lookup(
    haplo_frame: DataFrame, transcript: str | None = None
) -> dict[str, dict[str, str]]:
    """
    build_haplosaurus_lookup builds a lookup dictionary for Haplosaurus annotations
    keyed by sample.

    Parameters
    ----------
    haplo_frame : DataFrame
        DataFrame containing Haplosaurus annotations.
    transcript : str | None, optional
        Transcript ID to filter the Haplosaurus annotations, by default None.

    Returns
    -------
    dict[str, dict[str, str]]
        Dictionary keyed by sample, with values being dictionaries of Haplosaurus
        annotations.

    Raises
    ------
    ValueError
        Raised if a sample occurs in multiple Haplosaurus records.
    """
    if transcript is not None:
        haplo_frame = haplo_frame.loc[
            haplo_frame["transcript"].eq(transcript)
        ].copy()  # type: ignore

    haplo_lookup: dict[
        str,
        dict[str, str],
    ] = {}

    for _, row in haplo_frame.iterrows():
        samples = _extract_haplo_samples(_clean_string(row["samples"]))

        annotation = {
            "haplo_hgvs_c": _clean_string(row["cds_haplotype"]),
            "haplo_hgvs_p": _clean_string(row["protein_haplotype"]),
            "haplo_consequence": _clean_string(row["cds_flags"]),
            "haplo_protein_flags": _clean_string(row["protein_flags"]),
            "haplo_variants": _clean_string(row["contributing_variants"]),
        }

        for sample in samples:
            if sample in haplo_lookup:
                raise ValueError(
                    f"Sample {sample!r} occurs in multiple " "Haplosaurus records."
                )

            haplo_lookup[sample] = annotation
    return haplo_lookup


def build_annotation_rows(
    sequence_frame: DataFrame,
    sequence_variants: dict[str, list[tuple[str, int, str, str]]],
    sample_column: str,
    vep_lookup: dict[str, dict[str, str]],
    haplo_lookup: dict[str, dict[str, str]],
    rename_dict: dict[str, str] | None = None,
    ambiguous_column: str = "n_ambigous_variants",
) -> list[dict[str, str]]:
    """
    Build annotation rows for each variant in the sequence frame.

    Parameters
    ----------
    sequence_frame : DataFrame
        DataFrame containing sequence information.
    sequence_variants : dict[str, list[tuple[str, int, str, str]]]
        Dictionary keyed by sample ID, with values being lists of tuples representing
        variants. Each tuple contains (chromosome, position, reference allele, alternate allele).
    sample_column : str
        Name of the column containing sample IDs.
    vep_lookup : dict[str, dict[str, str]]
        Dictionary keyed by sample ID, with values being dictionaries of VEP annotations.
    haplo_lookup : dict[str, dict[str, str]]
        Dictionary keyed by sample ID, with values being dictionaries of Haplosaurus annotations.

    Returns
    -------
    list[dict[str, str]]
        List of annotation rows for each variant.
    """
    # ------------------------------------------------------------
    # Create one annotation row per sequence
    # ------------------------------------------------------------

    annotation_rows: list[dict[str, Any]] = []

    for sample_value in sequence_frame[sample_column]:
        sample = str(sample_value)
        ambigous_n = sequence_frame.loc[
            sequence_frame[sample_column] == sample, ambiguous_column
        ].values[0]
        variants = sequence_variants.get(sample, [])
        n_variants = len(variants)

        # --------------------------------------------------------
        # Initialize annotation
        # --------------------------------------------------------

        annotation: dict[str, Any] = {
            sample_column: sample,
            "n_variants": n_variants,
            # All variants belonging to this sequence.
            "variants": ";".join(_variant_string(v) for v in variants),
            # ------------------------------------------------
            # VEP / HGVS
            # ------------------------------------------------
            "hgvs_g": "",
            "hgvs_c": "",
            "hgvs_p": "",
            "hgvs_g (all)": "",
            "hgvs_c (all)": "",
            "hgvs_p (all)": "",
            # ------------------------------------------------
            # Haplosaurus
            # ------------------------------------------------
            "haplo_hgvs_c": "",
            "haplo_hgvs_p": "",
            "haplo_consequence": "",
            "haplo_protein_flags": "",
            "haplo_variants": "",
            # ------------------------------------------------
            # Status
            # ------------------------------------------------
            "vep_status": "empty",
            "haplo_status": "empty",
            # ------------------------------------------------
            # Overall source
            # ------------------------------------------------
            "annotation_source": "",
        }

        # --------------------------------------------------------
        # VEP
        #
        # VEP is only used for single-variant sequences.
        # --------------------------------------------------------

        if n_variants == 1:
            variant = variants[0]
            hgvs = get_hgvc_for_variant(variant, vep_lookup, rename_dict)

            if hgvs is not None:
                annotation["hgvs_g"] = hgvs[0]
                annotation["hgvs_c"] = hgvs[1]
                annotation["hgvs_p"] = hgvs[2]
                annotation["hgvs_g (all)"] = (annotation["hgvs_g"],)
                annotation["hgvs_c (all)"] = (annotation["hgvs_c"],)
                annotation["hgvs_p (all)"] = (annotation["hgvs_p"],)

                # The record exists, even if some individual
                # HGVS fields happen to be empty.
                annotated = any(
                    (
                        annotation["hgvs_g"],
                        annotation["hgvs_c"],
                        annotation["hgvs_p"],
                    )
                )
                annotation["vep_status"] = "annotated" if annotated else "not_found"

        elif n_variants > 1:
            # VEP is deliberately not used as the primary
            # annotation for multi-variant haplotypes.
            #
            # Individual VEP consequences do not represent
            # the combined haplotype.
            collect_hgvsg = []
            collect_hgvsc = []
            collect_hgvsp = []
            for variant in variants:
                hgvs = get_hgvc_for_variant(variant, vep_lookup, rename_dict)
                if hgvs is not None:
                    if hgvs[0]:
                        collect_hgvsg.append(hgvs[0])
                    if hgvs[1]:
                        collect_hgvsc.append(hgvs[1])
                    if hgvs[2]:
                        collect_hgvsp.append(hgvs[2])
            annotation["hgvs_g (all)"] = ";".join(collect_hgvsg)
            annotation["hgvs_c (all)"] = ";".join(collect_hgvsc)
            annotation["hgvs_p (all)"] = ";".join(collect_hgvsp)
            annotation["vep_status"] = "not_applicable"

        elif n_variants == 0:
            # No variants present in the sequence.
            if ambigous_n > 0:
                status = "only ambigous variants"
            else:
                status = "reference-like"
            annotation["vep_status"] = status
            annotation["haplo_status"] = status
            annotation["haplo_consequence"] = "REF"

        # Haplosaurus is checked for EVERY sequence, including single-variant sequences
        haplo_annotation = haplo_lookup.get(sample)

        if haplo_annotation is not None:
            annotation["haplo_hgvs_c"] = rename_annotation(
                haplo_annotation["haplo_hgvs_c"],
                rename_dict=rename_dict,
            )

            annotation["haplo_hgvs_p"] = rename_annotation(
                haplo_annotation["haplo_hgvs_p"],
                rename_dict=rename_dict,
            )

            annotation["haplo_consequence"] = haplo_annotation["haplo_consequence"]

            annotation["haplo_protein_flags"] = haplo_annotation["haplo_protein_flags"]

            annotation["haplo_variants"] = rename_annotation(
                haplo_annotation["haplo_variants"],
                rename_dict=rename_dict,
            )

            annotated = any(
                (
                    annotation["haplo_hgvs_c"],
                    annotation["haplo_hgvs_p"],
                    annotation["haplo_consequence"],
                    annotation["haplo_variants"],
                )
            )
            annotation["haplo_status"] = "annotated" if annotated else "empty"

        else:
            annotation["haplo_status"] = "not_found"

        # --------------------------------------------------------
        # Determine combined annotation source
        # --------------------------------------------------------

        vep_annotated = annotation["vep_status"] == "annotated"

        haplo_annotated = annotation["haplo_status"] == "annotated"

        annotation["annotation_source"] = get_annotation_source(
            vep_annotated=vep_annotated,
            haplo_annotated=haplo_annotated,
        )

        annotation_rows.append(annotation)
    return annotation_rows


def get_hgvc_for_variant(
    variant, vep_lookup, rename_dict
) -> tuple[str, str, str] | None:
    uploaded_variation = _uploaded_variation_from_variant(variant)

    vep_annotation = vep_lookup.get(uploaded_variation)

    if vep_annotation is not None:
        hgvs_g = rename_annotation(
            vep_annotation["hgvs_g"],
            rename_dict=rename_dict,
        )

        hgvs_c = rename_annotation(
            vep_annotation["hgvs_c"],
            rename_dict=rename_dict,
        )

        hgvs_p = rename_annotation(
            vep_annotation["hgvs_p"],
            rename_dict=rename_dict,
        )
        return hgvs_g, hgvs_c, hgvs_p
    return None


def get_annotation_source(vep_annotated: bool, haplo_annotated: bool) -> str:
    if vep_annotated and haplo_annotated:
        return "VEP+Haplosaurus"

    elif vep_annotated:
        return "VEP"

    elif haplo_annotated:
        return "Haplosaurus"

    else:
        return "none"


def rename_annotation(value: Any, rename_dict: dict[str, str] | None = None) -> str:
    """
    Replace identifiers in an annotation string using rename_dict.

    Replacement is performed longest-first to avoid partial
    replacements when one identifier is a substring of another.
    """
    value = _clean_string(value)

    if not value or not rename_dict:
        return value

    result = value

    for old_id, new_id in sorted(
        rename_dict.items(),
        key=lambda item: len(str(item[0])),
        reverse=True,
    ):
        old_id = str(old_id)
        new_id = str(new_id)

        if old_id:
            result = result.replace(old_id, new_id)

    return result

    # @element
    # def build_input(
    #     self,
    #     variants: Element | FileSource,
    #     *,
    #     wt_sequence: str | None = None,
    #     input_format: str = "ensembl",
    #     default_chrom: str = "1",
    #     default_genomic_start: int = 1,
    #     default_strand: str = "+",
    #     tag: TagType | None = None,
    #     out: OutType | None = None,
    #     par: ParType | None = None,
    #     cfg: CfgType | None = None,
    # ) -> Element:
    #     """Build a VEP input file from a variant/sequence table.

    #     The input table may either:
    #     1) already contain variant-level columns (``genomic_pos``, ``ref``,
    #        ``alt``), or
    #     2) contain sequence rows (mutated sequence + WT sequence).

    #     Column names can be overridden via ``par`` keys:
    #     ``sequence_column``, ``wt_column``, ``chrom_column``,
    #     ``start_column``, ``strand_column``, ``id_column``, ``ref_column``,
    #     ``alt_column``, ``pos_column``.
    #     """
    #     _ = cfg  # kept for API symmetry
    #     par = par or Params()
    #     fmt = input_format.lower()
    #     ext = "vcf" if fmt == "vcf" else "txt"

    #     tag = variants.tag.bump(
    #         stage=Stage.PREP,
    #         method=Method.ENSEMBL,
    #         state=State.GENERATED,
    #     ).resolve(tag)
    #     out = OutSpec(
    #         stem=f"{variants.tag.root}.vep_input.{fmt}",
    #         folder=self.default_outdir(variants.tag.root) / "input",
    #         ext=ext,
    #     ).resolve(out)

    #     runner = self.build_input_file(
    #         source_file=variants.file,
    #         output_file=out.file,
    #         wt_sequence=wt_sequence,
    #         input_format=fmt,
    #         default_chrom=default_chrom,
    #         default_genomic_start=default_genomic_start,
    #         default_strand=default_strand,
    #         par=par,
    #     )
    #     key = Element.generate_key(tag, "build_input", input_format=fmt)
    #     determinants = (
    #         f"input_format={fmt}",
    #         f"default_chrom={default_chrom}",
    #         f"default_genomic_start={default_genomic_start}",
    #         f"default_strand={default_strand}",
    #         f"has_wt_sequence={wt_sequence is not None}",
    #     ) + par.determinants

    #     return Element(
    #         key,
    #         runner,
    #         tag=tag,
    #         artifacts=ArtifactSet.from_outspec(out),
    #         determinants=determinants,
    #         inputs=(variants.file,),
    #         pres=(variants,),
    #     )

    # def build_input_file(
    #     self,
    #     source_file: Path,
    #     output_file: Path,
    #     *,
    #     wt_sequence: str | None,
    #     input_format: str,
    #     default_chrom: str,
    #     default_genomic_start: int,
    #     default_strand: str,
    #     par: ParType,
    # ) -> Runnable:
    #     """Create a runnable that writes a VEP input file."""

    #     def run() -> None:
    #         df = read_frame(source_file)

    #         sequence_col = str(par.get("sequence_column", "sequence"))
    #         wt_col = str(par.get("wt_column", "wt_sequence"))
    #         chrom_col = str(par.get("chrom_column", "chromosome"))
    #         start_col = str(par.get("start_column", "genomic_start"))
    #         strand_col = str(par.get("strand_column", "strand"))
    #         id_col = str(par.get("id_column", "seq_id"))
    #         ref_col = str(par.get("ref_column", "ref"))
    #         alt_col = str(par.get("alt_column", "alt"))
    #         pos_col = str(par.get("pos_column", "genomic_pos"))

    #         lines: list[str] = []
    #         for _, row in df.iterrows():
    #             chrom = _clean_string(row.get(chrom_col), default=default_chrom)
    #             strand = _to_strand_token(
    #                 _clean_string(row.get(strand_col), default=default_strand)
    #             )
    #             variant_id = None
    #             if id_col in row and not _is_missing(row.get(id_col)):
    #                 variant_id = str(row.get(id_col))

    #             row_variant: _SimpleVariant | None = None
    #             row_wt_seq: str | None = None
    #             row_genomic_start: int | None = None

    #             has_direct = (
    #                 pos_col in row
    #                 and ref_col in row
    #                 and alt_col in row
    #                 and not _is_missing(row.get(pos_col))
    #                 and not _is_missing(row.get(ref_col))
    #                 and not _is_missing(row.get(alt_col))
    #             )
    #             if has_direct:
    #                 start = int(row.get(pos_col))
    #                 ref = str(row.get(ref_col)).upper()
    #                 alt = str(row.get(alt_col)).upper()
    #                 end = start + len(ref) - 1
    #                 if ref == "-":
    #                     end = start - 1
    #                 row_variant = _SimpleVariant(
    #                     chrom=chrom,
    #                     start=start,
    #                     end=end,
    #                     ref=ref,
    #                     alt=alt,
    #                     strand=strand,
    #                     variant_id=variant_id,
    #                 )
    #                 row_genomic_start = int(row.get(start_col, default_genomic_start))
    #                 if wt_col in row and not _is_missing(row.get(wt_col)):
    #                     row_wt_seq = str(row.get(wt_col)).upper()
    #                 elif wt_sequence is not None:
    #                     row_wt_seq = wt_sequence.upper()
    #             else:
    #                 if sequence_col not in row or _is_missing(row.get(sequence_col)):
    #                     continue

    #                 mutated = str(row.get(sequence_col)).upper()
    #                 if wt_col in row and not _is_missing(row.get(wt_col)):
    #                     row_wt_seq = str(row.get(wt_col)).upper()
    #                 elif wt_sequence is not None:
    #                     row_wt_seq = wt_sequence.upper()

    #                 if row_wt_seq is None:
    #                     raise ValueError(
    #                         "No WT sequence available. Provide wt_sequence or a WT column in the input table."  # noqa: E501
    #                     )

    #                 row_genomic_start = int(row.get(start_col, default_genomic_start))
    #                 row_variant = _diff_variant_from_sequences(
    #                     row_wt_seq,
    #                     mutated,
    #                     row_genomic_start,
    #                 )
    #                 if row_variant is None:
    #                     continue
    #                 row_variant = _SimpleVariant(
    #                     chrom=chrom,
    #                     start=row_variant.start,
    #                     end=row_variant.end,
    #                     ref=row_variant.ref,
    #                     alt=row_variant.alt,
    #                     strand=strand,
    #                     variant_id=variant_id,
    #                 )

    #             line = _render_line(
    #                 row_variant,
    #                 input_format,
    #                 wt_sequence=row_wt_seq,
    #                 genomic_start=row_genomic_start,
    #             )
    #             lines.append(line)

    #         output_file.parent.mkdir(parents=True, exist_ok=True)
    #         with output_file.open("w", encoding="utf-8") as handle:
    #             if input_format == "vcf":
    #                 handle.write("##fileformat=VCFv4.2\n")
    #                 handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    #             if lines:
    #                 handle.write("\n".join(lines) + "\n")

    #     return Runnable(
    #         run,
    #         display=(
    #             f"vep.build_input_file(source={source_file}, out={output_file}, "
    #             f"format={input_format})"
    #         ),
    #     )


class Haplo(External):
    """Wrapper for the ``haplo`` CLI."""

    def __init__(
        self,
        name: str = "haplo",
        primary_binary: str = "haplo",
        version: str | None = None,
        source: str = "https://www.ensembl.org/info/docs/tools/vep/",
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
    ) -> None:
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=parameters or {},
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        """Detect HAPLO version from ``haplo --help`` output."""
        if self._version:
            return self._version
        if not self.primary_binary or not self.ensure_binary():
            return fallback
        try:
            cp = subprocess.run(
                [self.primary_binary, "--help"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            out = (cp.stdout or "") + "\n" + (cp.stderr or "")
            m = re.search(r"haplo\s*:\s*([0-9]+(?:\.[0-9]+)*)", out)
            if m:
                return m.group(1)
        except Exception:
            pass
        return fallback

    def default_outdir(self, sample_name: str) -> Path:
        """Default result location for VEP output files."""
        return Path("results") / "vep" / Vep().version_name / sample_name

    @property
    def haplo_params(self) -> Params:
        return Params()

    @element
    def haplo(
        self,
        haplo_input: Element | FileSource,
        *,
        input_format: str = "vcf",
        species: str = "mouse",
        transcript_id: str | None = None,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Run VEP for one prepared input file."""
        params = self.haplo_params.resolve(par)  # type: ignore
        tag = haplo_input.tag.bump(
            stage=Stage.ANALYSIS,
            method=Method.HAPLO,
            state=State.ANNOTATED,
            flag=species,
        ).resolve(tag)
        out = OutSpec.from_tag(
            tag, Out(folder=self.default_outdir(tag.root), ext="txt")
        ).resolve(out)
        artifacts = ArtifactSet.from_outspec(out)
        runner = self.run_haplo(
            input_file=haplo_input.file,
            output_file=out.file,
            input_format=input_format,
            transcript_id=transcript_id,
            par=par,
            cfg=cfg,
        )
        key = Element.generate_key(
            tag,
            "run",
            "haplo",
        )
        determinants = (species, input_format) + params.determinants

        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=(haplo_input.file,),
            pres=(haplo_input,),
            empty_ok=True,
        )

    @subroutine
    def run_haplo(
        self,
        input_file: Path,
        output_file: Path,
        *,
        input_format: str = "vcf",
        species: str = "mouse",
        transcript_id: str | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for one ``haplo`` invocation."""
        _ = cfg  # consumed by @subroutine wrapper
        par = self.haplo_params.resolve(par)  # type: ignore

        def post(output_file=output_file):
            if not output_file.exists():
                output_file.touch()

        vep_format = "ensembl" if input_format == "default" else input_format
        arguments = [
            "--input_file",
            str(Path(input_file).resolve()),
            "--output_file",
            str(Path(output_file).resolve()),
            "--species",
            str(species),
            "--format",
            str(vep_format),
            "--force_overwrite",
            "--database",
        ]
        if transcript_id is not None:
            arguments.extend(["--transcript", str(transcript_id)])
        arguments.extend(
            _params_to_cli(
                par,
            )
        )

        in_paths = [Path(input_file).resolve()]
        out_paths = [Path(output_file).resolve()]
        return (
            arguments,
            None,
            in_paths,
            out_paths,
            None,
            None,
            post,
        )


def variant_type_from_alleles(
    ref: str,
    alt: str,
) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"

    if len(ref) == len(alt):
        return "MNV"

    if len(ref) < len(alt):
        return "INS"

    if len(ref) > len(alt):
        return "DEL"

    raise ValueError(f"Could not determine variant type: {ref}>{alt}")

    @element
    def cache(self) -> Element:
        """build an offline cache for haplo"""
        pass

    @element
    def haplo(self) -> Element:
        """run haplo on the cache"""
        pass


# workflow
# for each sequence
# map to reference (single reference)
# turn alignment to edits
# normalize edits?
# create a list of Variants
#
#


# def create_vcf(
#     sequences: list[SequenceRecord],
#     reference: SequenceRecord,
#     chrom: str,
#     reference_start: int,
#     output_vcf: str | Path,
#     aligner: PairwiseAligner,
#     ploidy: int = 1,
#     allow_ambiguous: bool = False,
#     sample_name_prefix: str = "",
# ) -> dict:
#     """
#     Create a multi-sample VCF from consensus sequences.

#     Each FASTA record becomes one VCF sample.

#     For haploid output:
#         reference = (0,)
#         alternate = (1,)
#         missing   = (None,)

#     For diploid output:
#         reference = (0, 0)
#         alternate = (1, 1)
#         missing   = (None, None)

#     Note:
#         Consensus sequences cannot tell us whether a diploid biological
#         sample is heterozygous. Therefore ploidy=2 should only be used
#         if the biological interpretation justifies treating the consensus
#         as homozygous.
#     """
#     if ploidy not in (1, 2):
#         raise ValueError("Only ploidy=1 and ploidy=2 are supported")

#     output_vcf = Path(output_vcf)

#     if output_vcf.suffix != ".gz":
#         raise ValueError("Output VCF should end with .vcf.gz")

#     reference_length = len(reference.sequence)

#     # ------------------------------------------------------------------
#     # Call variants for all sequences.
#     # ------------------------------------------------------------------
#     sample_variants: dict[str, set[tuple]] = {}

#     stats = Counter()
#     sequence_stats = {}

#     used_sample_names: set[str] = set()

#     for sequence_record in sequences:
#         sample_name = sample_name_prefix + sequence_record.id

#         if sample_name in used_sample_names:
#             raise ValueError(f"Duplicate VCF sample name: {sample_name}")

#         used_sample_names.add(sample_name)

#         variants, seq_stats = call_variants_for_sequence(
#             reference=reference.sequence,
#             query=sequence_record.sequence,
#             chrom=chrom,
#             reference_start=reference_start,
#             aligner=aligner,
#             allow_ambiguous=allow_ambiguous,
#         )

#         variant_keys = set()

#         for variant in variants:
#             key = (
#                 variant.chrom,
#                 variant.pos,
#                 variant.ref,
#                 variant.alt,
#             )

#             variant_keys.add(key)
#             stats[variant.variant_type] += 1

#         sample_variants[sample_name] = variant_keys
#         sequence_stats[sample_name] = seq_stats

#         stats["sequences"] += 1

#         if not variants:
#             stats["identical_to_reference"] += 1
#         else:
#             stats["sequences_with_variants"] += 1

#     # ------------------------------------------------------------------
#     # Aggregate all variants.
#     # ------------------------------------------------------------------
#     all_variant_keys = set()

#     for variant_keys in sample_variants.values():
#         all_variant_keys.update(variant_keys)

#     sorted_variants = sorted(
#         all_variant_keys,
#         key=lambda x: (
#             x[0],
#             x[1],
#             x[2],
#             x[3],
#         ),
#     )

#     # ------------------------------------------------------------------
#     # VCF header.
#     # ------------------------------------------------------------------
#     header = pysam.VariantHeader()

#     header.add_meta(
#         "fileformat",
#         value="VCFv4.3",
#     )

#     header.add_meta(
#         "source",
#         value="amplicon_to_vcf",
#     )

#     header.add_meta(
#         "reference",
#         value=reference.id,
#     )

#     header.add_meta(
#         "contig",
#         items=[
#             ("ID", chrom),
#             ("length", reference_length),
#         ],
#     )

#     header.add_meta(
#         "INFO",
#         items=[
#             ("TYPE", "String", "1", "Variant type: SNV/MNV/INS/DEL"),
#         ],
#     )

#     header.add_meta(
#         "FORMAT",
#         items=[
#             ("GT", "String", "1", "Genotype"),
#         ],
#     )

#     for sequence_record in sequences:
#         sample_name = sample_name_prefix + sequence_record.id
#         header.add_sample(sample_name)

#     # ------------------------------------------------------------------
#     # Write compressed VCF.
#     # ------------------------------------------------------------------
#     with pysam.VariantFile(
#         str(output_vcf),
#         mode="wz",
#         header=header,
#     ) as vcf:
#         for chrom_, pos, ref, alt in sorted_variants:
#             variant = Variant(
#                 chrom=chrom_,
#                 pos=pos,
#                 ref=ref,
#                 alt=alt,
#                 variant_type=variant_type_from_alleles(ref, alt),
#             )

#             record = vcf.new_record(
#                 contig=chrom_,
#                 start=pos - 1,
#                 stop=pos - 1 + len(ref),
#                 alleles=(ref, alt),
#             )

#             record.info["TYPE"] = variant.variant_type

#             key = (
#                 variant.chrom,
#                 variant.pos,
#                 variant.ref,
#                 variant.alt,
#             )

#             for sequence_record in sequences:
#                 sample_name = sample_name_prefix + sequence_record.id

#                 if key in sample_variants[sample_name]:
#                     gt = (1,) if ploidy == 1 else (1, 1)
#                 else:
#                     gt = (0,) if ploidy == 1 else (0, 0)

#                 record.samples[sample_name]["GT"] = gt

#             vcf.write(record)

#     # ------------------------------------------------------------------
#     # Create tabix index.
#     # ------------------------------------------------------------------
#     pysam.tabix_index(
#         str(output_vcf),
#         preset="vcf",
#         force=True,
#     )

#     return {
#         "n_sequences": stats["sequences"],
#         "n_sequences_with_variants": stats["sequences_with_variants"],
#         "n_identical_to_reference": stats["identical_to_reference"],
#         "n_unique_variants": len(sorted_variants),
#         "snv_count": stats["SNV"],
#         "mnv_count": stats["MNV"],
#         "insertion_count": stats["INS"],
#         "deletion_count": stats["DEL"],
#         "sequence_stats": sequence_stats,
#     }


################################################################################
# combine haplo and vep
################################################################################


def _parse_vep_extra(extra: Any) -> dict[str, str]:
    """Parse the semicolon-separated VEP Extra column."""
    if pd.isna(extra):
        return {}

    result: dict[str, str] = {}

    for item in str(extra).split(";"):
        if "=" not in item:
            continue

        key, value = item.split("=", 1)
        result[key] = value

    return result


def _read_vep_output(path: Path) -> pd.DataFrame:
    """Read default tab-delimited VEP output."""
    with open(path) as f:
        header = next(line for line in f if line.startswith("#Uploaded_variation"))

    columns = header.rstrip("\n").lstrip("#").split("\t")

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        names=columns,
    )

    return df


def _read_haplo_output(path: Path) -> pd.DataFrame:
    """Read default 8-column Haplosaurus output."""
    columns = [
        "transcript",
        "cds_haplotype",
        "cds_flags",
        "protein_haplotype",
        "protein_flags",
        "frequency",
        "contributing_variants",
        "samples",
    ]

    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=columns,
        comment="#",
        dtype=str,
        keep_default_na=False,
    )


def _extract_haplo_samples(value: str) -> list[str]:
    """
    Extract sample IDs from the Haplosaurus sample:count field.

    Example:
        'seq_001:1,seq_002:1'
        -> ['seq_001', 'seq_002']
    """
    if not value:
        return []

    samples: list[str] = []

    for item in value.split(","):
        item = item.strip()

        if not item:
            continue

        # Haplosaurus uses sample:count.
        # split only once because sample IDs may contain ':'.
        sample, _, _count = item.rpartition(":")

        if not sample:
            # Defensive fallback if no count is present.
            sample = item

        samples.append(sample)

    return samples


def _select_vep_transcript(
    frame: DataFrame,
    transcript: str | None,
) -> DataFrame:
    """Restrict VEP annotations to the requested transcript."""
    if transcript is None:
        return frame

    return cast(DataFrame, frame.loc[frame["Feature"] == transcript, :])


def _build_vep_variant_annotations(
    vep_frame: pd.DataFrame,
    *,
    transcript: str | None,
) -> dict[str, dict[str, str]]:
    """
    Build lookup:

        uploaded_variant -> {hgvs_g, hgvs_c, hgvs_p}

    VEP default output contains one row per transcript/consequence.
    """
    vep_frame = _select_vep_transcript(
        vep_frame,
        transcript=transcript,
    )

    annotations: dict[str, dict[str, str]] = {}

    for _, row in vep_frame.iterrows():
        uploaded_variant = str(row["Uploaded_variation"])

        extra = _parse_vep_extra(row.get("Extra"))

        annotation = {
            "hgvs_g": extra.get("HGVSg", ""),
            "hgvs_c": extra.get("HGVSc", ""),
            "hgvs_p": extra.get("HGVSp", ""),
        }

        # Keep the first annotation for a given variant.
        # With a transcript filter there should normally be only one.
        annotations.setdefault(
            uploaded_variant,
            annotation,
        )

    return annotations


def _build_single_variant_annotations(
    sequence_frame: pd.DataFrame,
    vep_annotations: dict[str, dict[str, str]],
    *,
    sample_column: str,
    chrom_column: str,
    position_column: str,
    ref_column: str,
    alt_column: str,
) -> dict[str, dict[str, str]]:
    """
    Build final VEP annotations for sequences carrying exactly one variant.

    The key is the original sequence/sample ID.
    """
    result: dict[str, dict[str, str]] = {}

    grouped = sequence_frame.groupby(sample_column, sort=False)

    for sample, group in grouped:
        sample = str(sample)

        if len(group) != 1:
            continue

        row = group.iloc[0]

        chrom = str(row[chrom_column])
        position = int(row[position_column])
        ref = str(row[ref_column]).upper()
        alt = str(row[alt_column]).upper()

        # This must match VEP's Uploaded_variation representation.
        uploaded_variant = f"{chrom}_{position}_{ref}/{alt}"

        annotation = vep_annotations.get(uploaded_variant)

        if annotation is None:
            continue

        result[sample] = annotation

    return result


def _build_haplo_annotations(
    haplo_frame: DataFrame,
    *,
    transcript: str | None,
) -> dict[str, dict[str, str]]:
    """
    Map Haplosaurus haplotypes back to individual samples.

    Returns:

        sample -> {
            haplo_hgvs_c: ...,
            haplo_hgvs_p: ...,
            haplo_consequence: ...,
            haplo_flags_c: ...,
            haplo_flags_p: ...,
            haplo_variants: ...
        }
    """
    if transcript is not None:
        haplo_frame = cast(
            DataFrame,
            haplo_frame.loc[
                haplo_frame["transcript"].astype(str).eq(transcript)
            ].copy(),
        )

    result: dict[str, dict[str, str]] = {}

    for _, row in haplo_frame.iterrows():
        samples = _extract_haplo_samples(str(row["samples"]))

        annotation = {
            "haplo_hgvs_c": str(row["cds_haplotype"]),
            "haplo_hgvs_p": str(row["protein_haplotype"]),
            "haplo_consequence": str(row["cds_flags"]),
            "haplo_protein_flags": str(row["protein_flags"]),
            "haplo_variants": str(row["contributing_variants"]),
        }

        for sample in samples:
            # A sample can theoretically occur in multiple haplotype
            # rows, so fail loudly rather than silently overwriting.
            if sample in result:
                raise ValueError(
                    f"Sample {sample!r} occurs in multiple Haplosaurus "
                    "haplotype records for the selected transcript."
                )

            result[sample] = annotation.copy()

    return result


def _variant_key_from_row(
    row: Series,
    *,
    chrom_column: str,
    position_column: str,
    ref_column: str,
    alt_column: str,
) -> tuple[str, int, str, str]:
    """Create a normalized variant key."""
    return (
        _clean_string(row[chrom_column]),
        int(row[position_column]),
        _clean_string(row[ref_column]).upper(),
        _clean_string(row[alt_column]).upper(),
    )


def _variant_string(
    variant: tuple[str, int, str, str],
) -> str:
    """Convert a variant tuple into a compact human-readable string."""
    chrom, position, ref, alt = variant
    return f"{chrom}:{position}{ref}>{alt}"


def _uploaded_variation_from_variant(
    variant: tuple[str, int, str, str],
) -> str:
    """Create the VEP Uploaded_variation identifier."""
    chrom, position, ref, alt = variant
    return f"{chrom}_{position}_{ref}/{alt}"
