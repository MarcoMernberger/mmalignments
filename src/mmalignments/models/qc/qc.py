"""Convenience functions for post-mapping quality control workflows.

This module provides high-level functions that chain multiple QC steps
into a single workflow, similar to the alignsort pattern in BWAMem2.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Mapping

from mmalignments.models.aligners.samtools import Samtools
from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.callers.gatk import GATK
from mmalignments.models.data import Genome
from mmalignments.models.elements import (
    Element,
    FileSource,
    MappedElement,
    NextGenSample,
)
from mmalignments.models.overlay import (
    CfgType,
    ExternalRunConfig,
    OutSpec,
    OutType,
    PartialOutSpec,
    ParType,
    SpecValues,
    TagType,
    pick_spec,
)
from mmalignments.models.qc.fastp import FastP
from mmalignments.models.qc.fastqc import FastQC
from mmalignments.models.qc.mosdepth import Mosdepth
from mmalignments.models.qc.multiqc import MultiQC
from mmalignments.models.qc.parsers import build_qc_summary
from mmalignments.models.resources import ResourceConfig  # type: ignore[import]
from mmalignments.models.tags import Method, Stage, State
from mmalignments.services.io import (
    open_fastq,
    write_fastq_check_results,
    write_json,
)

from .parsers import parse_fastq_record

logger = logging.getLogger(__name__)

###############################################################################
# Configs
###############################################################################


@dataclass(frozen=True)
class PreQCConfig:
    qc_root: Path | None = None
    threads: int = field(default_factory=lambda: ResourceConfig.detect().threads)
    n_records: int = 10_000
    run_sanity_check: bool = True
    run_fastp: bool = True
    run_fastqc_raw: bool = True
    run_fastqc_cleaned: bool = True
    run_summary: bool = True


@dataclass(frozen=True)
class PostQCConfig:
    qc_root: Path | None = None
    threads: int = field(default_factory=lambda: ResourceConfig.detect().threads)
    run_alignment_summary: bool = True
    run_insert_size: bool = True
    run_hs_metrics: bool = True
    run_mosdepth: bool = True
    run_samtools_flagstat: bool = True
    run_samtools_stats: bool = True
    multiqc: bool = True


########################################################################################
# Helpers
########################################################################################


def _extract_read_id(header: str) -> str:
    """Extract read ID (first token, without leading '@', without /1 or /2)."""
    return header[1:].split()[0].split("/")[0]


def _validate_record(
    read_label: str,
    record_index: int,
    header: str,
    seq: str,
    plus: str,
    qual: str,
) -> list[str]:
    """
    Validate a single FASTQ record and return error messages (if any).

    Parameters
    ----------
    read_label : str
        Label for the read (e.g., "R1" or "R2")
    record_index : int
        Index of the record in the FASTQ file
    header : str
        Header line of the FASTQ record
    seq : str
        Sequence line of the FASTQ record
    plus : str
        Plus line of the FASTQ record
    qual : str
        Quality line of the FASTQ record

    Returns
    -------
    list[str]
        List of error messages, empty if no errors
    """
    errs: list[str] = []
    n = record_index + 1  # human-friendly

    if not header.startswith("@"):
        errs.append(f"{read_label} record {n}: Header doesn't start with '@'")
    if not plus.startswith("+"):
        errs.append(f"{read_label} record {n}: Plus line doesn't start with '+'")
    if len(seq) != len(qual):
        errs.append(
            f"{read_label} record {n}: Sequence and quality lengths don't match"
        )

    return errs


def _scan_fastq(
    path: Path | str,
    *,
    read_label: str,
    n_records: int,
) -> tuple[list[str], int, list[str]]:
    """
    Scan up to n_records from a FASTQ and return IDs, record count, and errors.

    Parameters
    ----------
    path : Path | str
        Path to the FASTQ file
    read_label : str
        Label for the read (e.g., "R1" or "R2")
    n_records : int
        Maximum number of records to scan

    Returns
    -------
    tuple[list[str], int, list[str]]
        Tuple containing:
        - List of read IDs
        - Number of records checked
        - List of error messages
    """
    ids: list[str] = []
    errors: list[str] = []
    count = 0

    with open_fastq(path) as fh:
        for i in range(n_records):
            record = parse_fastq_record(fh)
            if record is None:
                break

            header, seq, plus, qual = record
            count += 1

            errors.extend(_validate_record(read_label, i, header, seq, plus, qual))
            ids.append(_extract_read_id(header))

    return (ids, count, errors)


def _pairing_mismatches(
    r1_ids: list[str], r2_ids: list[str], *, max_examples: int = 5
) -> tuple[int, list[str]]:
    """
    Compare paired IDs and return mismatch count and warning examples.

    Parameters
    ----------
    r1_ids : list[str]
        List of read IDs from R1
    r2_ids : list[str]
        List of read IDs from R2
    max_examples : int, optional
        Maximum number of mismatch examples to include in warnings, by default 5

    Returns
    -------
    tuple[int, list[str]]
        Tuple containing:
        - Number of mismatches
        - List of warning messages
    """
    mismatches = 0
    warnings: list[str] = []

    for i, (id1, id2) in enumerate(zip(r1_ids, r2_ids)):
        if id1 != id2:
            mismatches += 1
            if mismatches <= max_examples:
                warnings.append(f"ID mismatch at record {i+1}: R1='{id1}' R2='{id2}'")

    return mismatches, warnings


def _validate_targets_baits(
    cfg: PostQCConfig,
    targets: Element | FileSource | None,
    baits: Element | FileSource | None,
) -> tuple[Element | FileSource | None, Element | FileSource | None]:
    """
    Validate and prepare targets and baits for HS metrics and mosdepth.

    Parameters
    ----------
    cfg : PostQCConfig
        Configuration for post-alignment QC.
    targets : Element | FileSource | None
        Target intervals for HS metrics and mosdepth.
    baits : Element | FileSource | None
        Bait intervals for HS metrics.

    Returns
    -------
    tuple[Element | FileSource | None, Element | FileSource | None]
        Validated and prepared target and bait intervals.

    Raises
    ------
    ValueError
        If required targets or baits are not provided when HS metrics or mosdepth are
        enabled.
    """
    if cfg.run_hs_metrics:
        if targets is None:
            raise ValueError(
                "run_hs_metrics=True requires either "
                "(bait_intervals and target_intervals) or targets to be provided"
            )
        # if baits is None:
        # logger.warning(
        #     "bait_intervals and/or target_intervals not provided, "
        #     "using targets for both. This may not be ideal for HS metrics."
        # )
        # baits = targets

    if cfg.run_mosdepth and targets is None:
        raise ValueError("run_mosdepth=True requires targets to be provided")

    return targets, baits


########################################################################################
# Elements
########################################################################################


def _step_out(
    out: OutType | SpecValues[OutType] | None, step: str, folder: Path | str
) -> OutType:
    """Return the output spec for ``step`` with ``folder`` as default folder."""
    return PartialOutSpec(folder=Path(folder)).resolve(pick_spec(out, step))


def build_input_check(
    sample: NextGenSample,
    n_records: int = 10000,
    *,
    tag: TagType | None = None,
    out: OutType | None = None,
) -> Element:
    """Create Element that performs FASTQ sanity checks."""
    resolved_tag = sample.tag.bump(
        stage=Stage.QC,
        method=Method.CUSTOM,
        state=State.STAT,
    ).resolve(tag)
    resolved_out = OutSpec.from_tag(resolved_tag).resolve(
        OutSpec(folder=sample.r1.parent / "qc" / "input_check", ext="json"),
        out,
    )

    out_json = resolved_out.file
    out_txt = out_json.with_suffix(".txt")

    def _check_fastq() -> None:
        results: dict[str, Any] = {
            "sample_name": sample.root,
            "pairing": sample.pairing.value,
            "checks": {},
            "errors": [],
            "warnings": [],
        }

        try:
            # R1 scan
            r1_ids, r1_counts, r1_errors = _scan_fastq(
                sample.r1, read_label="R1", n_records=n_records
            )
            results["checks"]["r1_records_checked"] = r1_counts
            results["errors"].extend(r1_errors)
            results["checks"]["r1_format_ok"] = len(r1_errors) == 0

            # R2 scan + pairing check (only if paired)
            if sample.pairing == "paired" and sample.r2:
                r2_ids, r2_counts, r2_errors = _scan_fastq(
                    sample.r2, read_label="R2", n_records=n_records
                )
                results["checks"]["r2_records_checked"] = r2_counts
                results["errors"].extend(r2_errors)
                results["checks"]["r2_format_ok"] = len(r2_errors) == 0

                mismatches, mismatch_warnings = _pairing_mismatches(r1_ids, r2_ids)
                results["checks"]["id_mismatches"] = mismatches
                results["checks"]["ids_match"] = mismatches == 0
                results["warnings"].extend(mismatch_warnings)

                if mismatches > 0:
                    results["warnings"].append(
                        f"Found {mismatches} ID mismatches in {min(r1_counts, r2_counts)} checked records"  # noqa: E501
                    )

            results["status"] = "PASS" if not results["errors"] else "FAIL"

        except Exception as e:
            results["status"] = "ERROR"
            results["errors"].append(f"Exception during check: {e}")
            logger.error(f"Error checking FASTQ files for {sample.root}: {e}")

        write_json(out_json, results)
        write_fastq_check_results(out_txt, sample.root, sample.pairing, results)

        logger.info(f"Input check completed for {sample.root}: {results['status']}")

    input_files = [sample.r1] + ([sample.r2] if sample.r2 else [])
    return Element(
        Element.generate_key(
            resolved_tag, "custom", "build_input_check", n_records=n_records
        ),
        _check_fastq,
        tag=resolved_tag,
        determinants=(f"n_records={n_records}",),
        inputs=tuple(input_files),
        artifacts=ArtifactSet(
            out_json, primary_name="input_check_json", input_check_txt=out_txt
        ),
        pres=(),
    )


###############################################################################
# PreAlignment QC workflow
###############################################################################


def pre_alignment_qc(
    sample: NextGenSample,
    *,
    tag: TagType | SpecValues[TagType] | None = None,
    out: OutType | SpecValues[OutType] | None = None,
    par: ParType | SpecValues[ParType] | None = None,
    cfg: CfgType | SpecValues[CfgType] | None = None,
    preqc_cfg: PreQCConfig | None = None,
) -> tuple[Element, ...]:
    """Run comprehensive pre-alignment quality control workflow.

    Orchestrates a complete QC workflow for FASTQ files including:
    1. Input sanity checks (format validation, ID matching)
    2. fastp (quality filtering and trimming)
    3. FastQC on raw reads
    4. FastQC on cleaned reads
    5. QC summary (parse and aggregate metrics)

    Parameters
    ----------
    sample : NextGenSample
        Sample object containing FASTQ file paths and metadata.
    tag, out, par, cfg
        Overrides, optionally per step via ``SpecValues`` with the steps
        ``input_check``, ``fastp``, ``fastqc_raw``, ``fastqc_cleaned`` and
        ``summary``. Output folders default to
        ``preqc_cfg.qc_root`` or ``sample.result_dir / "qc"``.
    preqc_cfg : PreQCConfig | None
        Configuration object for pre-alignment QC. If None, defaults will
        be used.

    Returns
    -------
    tuple[Element, ...]
        Tuple of QC Elements in dependency order: input_check, fastp,
        fastqc_raw, fastqc_cleaned, qc_summary. The first element is the last
        in the dependency chain (qc_summary depends on all previous steps).

    Examples
    --------
    >>> from mmalignments.models.data import NextGenSample
    >>> from mmalignments.models.externals import ExternalRunConfig
    >>> sample = NextGenSample(root="sample1", pairing="paired",
    ...                        r1="raw_R1.fastq.gz",
    ...                        r2="raw_R2.fastq.gz")
    >>> qc_elements = pre_alignment_qc(
    ...     sample,
    ...     preqc_cfg=PreQCConfig(qc_root=Path("results/qc"), threads=8),
    ... )
    >>> for elem in qc_elements:
    ...     elem.run()
    """
    preqc_cfg = preqc_cfg or PreQCConfig()
    qc_dir = preqc_cfg.qc_root or (sample.result_dir / "qc")
    qc_elements = []
    default_cfg = ExternalRunConfig(threads=preqc_cfg.threads)

    # 1. Input check
    logger.info(f"Creating input check for {sample.root}")
    elements_to_summary: list[Element] = []
    if preqc_cfg.run_sanity_check:
        logger.info(f"Adding input sanity check for {sample.root}")
        input_check = build_input_check(
            sample=sample,
            n_records=preqc_cfg.n_records,
            tag=pick_spec(tag, "input_check"),
            out=_step_out(out, "input_check", qc_dir),
        )
        qc_elements.insert(0, input_check)

    # 2. fastp Element
    fastp_elem = None
    if preqc_cfg.run_fastp:
        logger.info(f"Creating fastp QC for {sample.root}")
        fastp = FastP()
        fastp_elem = fastp.qc(
            sample=sample,
            tag=pick_spec(tag, "fastp"),
            out=_step_out(out, "fastp", qc_dir),
            par=pick_spec(par, "fastp"),
            cfg=pick_spec(cfg, "fastp") or default_cfg,
        )
        qc_elements.insert(0, fastp_elem)

    # 3. FastQC on raw FASTQs
    fastqc = FastQC()
    if preqc_cfg.run_fastqc_raw:
        logger.info(f"Creating FastQC for raw reads of {sample.root}")
        fastqc_raw = fastqc.qc(
            sample=sample,
            label="raw",
            tag=pick_spec(tag, "fastqc_raw"),
            out=_step_out(out, "fastqc_raw", qc_dir),
            par=pick_spec(par, "fastqc_raw"),
            cfg=pick_spec(cfg, "fastqc_raw") or default_cfg,
        )
        qc_elements.insert(0, fastqc_raw)
        elements_to_summary.append(fastqc_raw)

    # 4. FastQC on cleaned FASTQs
    if preqc_cfg.run_fastqc_cleaned and fastp_elem:
        elements_to_summary.append(fastp_elem)
        logger.info(f"Creating FastQC for cleaned reads of {sample.root}")
        cleaned_sample = NextGenSample(
            root=f"{sample.root}_fastp",
            source=fastp_elem,
        )
        fastqc_cleaned = fastqc.qc(
            sample=cleaned_sample,  # Use fastp element as input (has cleaned FASTQs)
            label="cleaned",
            tag=pick_spec(tag, "fastqc_cleaned"),
            out=_step_out(out, "fastqc_cleaned", qc_dir),
            par=pick_spec(par, "fastqc_cleaned"),
            cfg=pick_spec(cfg, "fastqc_cleaned") or default_cfg,
        )
        qc_elements.insert(0, fastqc_cleaned)
        elements_to_summary.append(fastqc_cleaned)

    # 5. QC summary Element (pres=[fastp, fastqc_raw, fastqc_cleaned])
    if preqc_cfg.run_summary:
        logger.info(f"Creating QC summary for {sample.root}")
        qc_summary = build_qc_summary(
            sample=sample,
            elements_to_summary=tuple(elements_to_summary),
            tag=pick_spec(tag, "summary"),
            out=_step_out(out, "summary", qc_dir),
        )
        qc_elements.insert(0, qc_summary)

    logger.info(
        f"Created pre-alignment QC workflow with {len(qc_elements)} steps for {sample.root}"  # noqa: E501
    )
    return tuple(qc_elements)


########################################################################################
# PostAlignment QC workflow
########################################################################################


def post_mapping_qc_with_multiqc(
    mapped: MappedElement,
    reference: Genome,
    *,
    refdict_element: Element | FileSource | None = None,
    targets: Element | FileSource | None = None,
    baits: Element | FileSource | None = None,
    targets_intervals: Element | FileSource | None = None,
    bait_intervals: Element | FileSource | None = None,
    tag: TagType | SpecValues[TagType] | None = None,
    out: OutType | SpecValues[OutType] | None = None,
    par: ParType | SpecValues[ParType] | None = None,
    cfg: CfgType | SpecValues[CfgType] | None = None,
    main_qc_dir: Path | str | None = None,
    postqc_cfg: PostQCConfig | None = None,
) -> tuple[Element, Mapping[str, Element]]:
    """Run post-mapping QC and aggregate results with MultiQC.

    This is a convenience wrapper around post_mapping_qc() that also creates
    a MultiQC Element to aggregate all QC reports into a single HTML report.

    Parameters
    ----------
    mapped : MappedElement
        MappedElement containing the BAM file and metadata.
    reference : Genome
        Reference genome object with FASTA file path.
    targets : Element | FileSource | None
        BED file or Element containing target regions.
    baits : Element | FileSource | None
        Interval list for bait regions (Picard format).
    targets_intervals : Element | FileSource | None
        Interval list for target regions (Picard format).
    bait_intervals : Element | FileSource | None
        Interval list for bait regions (Picard format).
    tag, out, par, cfg
        Overrides, optionally per step via ``SpecValues``. Steps are those of
        :func:`post_mapping_qc` plus ``multiqc``. The ``multiqc`` output spec
        defines folder and file name of the aggregated HTML report.
        ``par`` for ``multiqc`` defaults to ``ParType(force=True)``.
    main_qc_dir : Path | str | None
        Base directory for QC outputs. If None, defaults to BAM parent / "qc".
    postqc_cfg : PostQCConfig | None, optional
        Configuration object for post-mapping QC, by default PostQCConfig()

    Returns
    -------
    tuple[Element, Mapping[str, Element]]
            Tuple of (multiqc_element, qc_elements) where qc_elements is the
            dictionary of individual QC Elements with keys indicating the step names.
            multiqc_element is the MultiQC aggregation Element that depends on all QC
            elements.
            The first returned element is the last Element in the QC chain and depends
            on the previous elements.

    Examples
    --------
    >>> multiqc_elem, qc_elems = post_mapping_qc_with_multiqc(
    ...     mapped=mapped_sample,
    ...     reference=genome,
    ...     targets=targets_bed,
    ...     baits=targets_interval_list,
    ... )
    >>> # Execute all QC steps
    >>> for elem in qc_elems.values():
    ...     elem.run()
    >>> # Generate MultiQC report
    >>> multiqc_elem.run()
    """
    postqc_cfg = postqc_cfg or PostQCConfig()
    if main_qc_dir is not None:
        postqc_cfg = replace(postqc_cfg, qc_root=Path(main_qc_dir))
    analysis_dir = postqc_cfg.qc_root or mapped.bam.parent / "qc"
    qc_elements = post_mapping_qc(
        mapped=mapped,
        reference=reference,
        refdict_element=refdict_element,
        targets=targets,
        baits=baits,
        targets_intervals=targets_intervals,
        bait_intervals=bait_intervals,
        tag=tag,
        out=out,
        par=par,
        cfg=cfg,
        postqc_cfg=postqc_cfg,
    )

    # Create MultiQC aggregation element
    multiqc = MultiQC()
    multiqc_elem = multiqc.aggregate(
        analysis_dir=analysis_dir,
        qc_elements=qc_elements,
        tag=pick_spec(tag, "multiqc"),
        out=pick_spec(out, "multiqc"),
        par=pick_spec(par, "multiqc") or ParType(force=True),
        cfg=pick_spec(cfg, "multiqc"),
    )

    logger.info(
        f"Created QC workflow with {len(qc_elements)} QC steps "
        f"and MultiQC aggregation for {mapped.name}"
    )

    return multiqc_elem, qc_elements


# post alignment qc
def post_mapping_qc(
    mapped: MappedElement,
    reference: Genome,
    *,
    refdict_element: Element | FileSource | None = None,
    targets: Element | FileSource | None = None,
    baits: Element | FileSource | None = None,
    targets_intervals: Element | FileSource | None = None,
    bait_intervals: Element | FileSource | None = None,
    tag: TagType | SpecValues[TagType] | None = None,
    out: OutType | SpecValues[OutType] | None = None,
    par: ParType | SpecValues[ParType] | None = None,
    cfg: CfgType | SpecValues[CfgType] | None = None,
    postqc_cfg: PostQCConfig | None = None,
) -> Mapping[str, Element]:
    """
    Run comprehensive post-mapping quality control on a mapped sample.

    This function orchestrates multiple QC tools to assess alignment quality,
    coverage, insert size distribution, and target enrichment. It creates a
    DAG of QC Elements that depend on the input mapped sample.

    The workflow includes:
    - Picard CollectAlignmentSummaryMetrics (mapping rates, mismatches)
    - Picard CollectInsertSizeMetrics (fragment length distribution)
    - Picard CollectHsMetrics (on-target rate, coverage uniformity)
    - mosdepth (depth distribution on targets)
    - samtools flagstat (basic alignment stats)
    - samtools stats (comprehensive alignment stats)

    Parameters
    ----------
    mapped : MappedElement
        MappedElement containing the BAM file and metadata.
    reference : Genome
        Reference genome object with FASTA file path.
    targets : Optional[Union[Element, Path, str]]
       BED file or Element containing target regions for mosdepth.
        Required if run_mosdepth=True or run_hs_metrics=True (and intervals not
        provided).
    baits : Optional[Union[Element, Path, str]]
        Interval list for bait regions (Picard format).
        Required if run_hs_metrics=True. If not provided but targets is,
        will attempt to use targets.
    targets_intervals : Optional[Union[Element, Path, str]]
        Interval list for target regions (Picard format).
    bait_intervals : Optional[Union[Element, Path, str]]
        Interval list for bait regions (Picard format).
        Run Picard CollectAlignmentSummaryMetrics (default: True).
    tag, out, par, cfg
        Overrides, optionally per step via ``SpecValues`` with the steps
        ``alignmetrics``, ``insertmetrics``, ``targets2interval``, ``hs``,
        ``mosdepth``, ``flagstat`` and ``stats``. Output folders default to
        ``picard``/``mosdepth``/``samtools`` below ``postqc_cfg.qc_root`` (or
        BAM parent / "qc").
    postqccfg : Optional[PostQCConfig], optional
        Configuration object for post-mapping QC, by default PostQCConfig()

    Returns
    -------
    tuple[Element, ...]
        Tuple of QC Elements that were created. Each Element represents
        a QC step and can be executed independently or as part of a DAG.

    Raises
    ------
    ValueError
        If required parameters are missing for enabled QC steps.

    Examples
    --------
    Run all QC steps on a mapped sample::

        >>> from mmalignments.models.qc.qc import post_mapping_qc
        >>> qc_elements = post_mapping_qc(
        ...     mapped=mapped_sample,
        ...     reference=genome,
        ...     targets=targets_bed,
        ...     bait_intervals=targets_interval_list,
        ...     target_intervals=targets_padded_interval_list,
        ...     threads_mosdepth=8
        ... )
        >>> for elem in qc_elements:
        ...     elem.run()

    Run only specific QC steps::

        >>> qc_elements = post_mapping_qc(
        ...     mapped=mapped_sample,
        ...     reference=genome,
        ...     run_alignment_summary=True,
        ...     run_insert_size=True,
        ...     run_hs_metrics=False,
        ...     run_mosdepth=False,
        ...     run_samtools_stats=True,
        ...     run_samtools_flagstat=True
        ... )

    Notes
    -----
    - All QC outputs are placed in subdirectories under the BAM file's parent
      directory by default (or output_base_dir if specified).
    - The function returns a tuple of Elements that form a DAG with the mapped
      sample as the root dependency.
    - For targeted sequencing (exomes, panels), provide interval lists for
      accurate on-target metrics.
    - mosdepth typically works with BED files, while Picard tools require
      interval_list format. Use GATK's BedToIntervalList to convert if needed.
    """
    postqc_cfg = postqc_cfg or PostQCConfig()
    output_dir = Path(postqc_cfg.qc_root or (mapped.bam.parent / "qc"))
    picard_dir = output_dir / "picard"
    mosdepth_dir = output_dir / "mosdepth"
    samtools_dir = output_dir / "samtools"

    qc_elements = {}

    # tools
    gatk = GATK()
    mosdepth = Mosdepth()
    st = Samtools()

    # Validate baits/targets for HS metrics and mosdepth
    targets, baits = _validate_targets_baits(postqc_cfg, targets, baits)

    # 1. Picard CollectAlignmentSummaryMetrics
    if postqc_cfg.run_alignment_summary:
        logger.info(f"Creating alignment summary metrics QC for {mapped.name}")
        alignment_elem = gatk.alignmetrics(
            mapped=mapped,
            reference=reference,
            tag=pick_spec(tag, "alignmetrics"),
            out=_step_out(out, "alignmetrics", picard_dir),
            par=pick_spec(par, "alignmetrics"),
            cfg=pick_spec(cfg, "alignmetrics"),
        )
        qc_elements["alignmetrics"] = alignment_elem

    # 2. Picard CollectInsertSizeMetrics
    if postqc_cfg.run_insert_size:
        logger.info(f"Creating insert size metrics QC for {mapped.name}")
        insert_elem = gatk.insertmetrics(
            mapped=mapped,
            tag=pick_spec(tag, "insertmetrics"),
            out=_step_out(out, "insertmetrics", picard_dir),
            par=pick_spec(par, "insertmetrics"),
            cfg=pick_spec(cfg, "insertmetrics"),
        )
        qc_elements["insertmetrics"] = insert_elem

    # 3. Picard CollectHsMetrics
    if postqc_cfg.run_hs_metrics:
        # create a sequence dict ... once!
        if refdict_element is None:
            raise ValueError(
                "run_hs_metrics=True requires a sequence dictionary for the reference. "
                "Please provide refdict_element or create one using GATK's CreateSequenceDictionary."  # noqa: E501
            )
        # create intervals from bed
        # targets_intervals = None
        # bait_intervals = None
        # if targets is not None:
        #     targets_intervals = gatk.bed2interval(
        #         bed_element=targets,
        #         sequence_dict=refdict_element,
        #         tag=pick_spec(tag, "targets2interval"),
        #         out=pick_spec(out, "targets2interval"),
        #         par=pick_spec(par, "targets2interval"),
        #         cfg=pick_spec(cfg, "targets2interval"),
        #     )
        # if bait_intervals is not None:
        # create another interval list for baits (Picard format)
        # bait_intervals = gatk.bed2interval(  # this must be outside
        #     bed_element=baits,
        #     sequence_dict=refdict_element,
        #     tag=pick_spec(tag, "targets2interval"),
        #     out=pick_spec(out, "targets2interval"),
        #     par=pick_spec(par, "targets2interval"),
        #     cfg=pick_spec(cfg, "targets2interval"),
        # )
        if bait_intervals is not None and targets_intervals is not None:
            logger.info(f"Creating HS metrics QC for {mapped.name}")
            hs_elem = gatk.hs(
                mapped=mapped,
                reference=reference,
                bait_intervals=bait_intervals,
                targets_intervals=targets_intervals,
                tag=pick_spec(tag, "hs"),
                out=_step_out(out, "hs", picard_dir),
                par=pick_spec(par, "hs"),
                cfg=pick_spec(cfg, "hs"),
            )
            qc_elements["hs"] = hs_elem

    # 4. mosdepth coverage
    if postqc_cfg.run_mosdepth:
        logger.info(f"Creating mosdepth coverage QC for {mapped.name}")
        # mosdepth.coverage adds the targets as ``by`` itself
        depth_elem = mosdepth.coverage(
            mapped=mapped,
            targets=targets,
            tag=pick_spec(tag, "mosdepth"),
            out=_step_out(out, "mosdepth", mosdepth_dir),
            par=pick_spec(par, "mosdepth") or ParType(threads=47, no_per_base=True),
            cfg=pick_spec(cfg, "mosdepth"),
        )
        qc_elements["mosdepth"] = depth_elem

    # 5. samtools flagstat
    if postqc_cfg.run_samtools_flagstat:
        logger.info(f"Creating samtools flagstat QC for {mapped.name}")
        flagstat_elem = st.flagstat(
            mapped=mapped,
            tag=pick_spec(tag, "flagstat"),
            out=_step_out(out, "flagstat", samtools_dir),
            par=pick_spec(par, "flagstat"),
            cfg=pick_spec(cfg, "flagstat"),
        )
        qc_elements["flagstat"] = flagstat_elem

    # 6. samtools stats
    if postqc_cfg.run_samtools_stats:
        logger.info(f"Creating samtools stats QC for {mapped.name}")
        stats_elem = st.stats(
            mapped=mapped,
            tag=pick_spec(tag, "stats"),
            out=_step_out(out, "stats", samtools_dir),
            par=pick_spec(par, "stats"),
            cfg=pick_spec(cfg, "stats"),
        )
        qc_elements["stats"] = stats_elem

    logger.info(f"Created {len(qc_elements)} QC elements for {mapped.name}")
    return qc_elements
