"""Convenience functions for post-mapping quality control workflows.

This module provides high-level functions that chain multiple QC steps
into a single workflow, similar to the alignsort pattern in BWAMem2.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

from mmalignments.models.aligners.samtools import Samtools
from mmalignments.models.callers.gatk import GATK
from mmalignments.models.data import Genome
from mmalignments.models.elements import Element, MappedElement, NextGenSampleElement
from mmalignments.models.externals import ExternalRunConfig
from mmalignments.models.parameters import Params
from mmalignments.models.qc.fastp import FastP
from mmalignments.models.qc.fastqc import FastQC
from mmalignments.models.qc.mosdepth import Mosdepth
from mmalignments.models.qc.multiqc import MultiQC
from mmalignments.models.qc.parsers import build_qc_summary
from mmalignments.models.tags import ElementTag, Method, Stage, State, merge_tag
from mmalignments.models.tags import PartialElementTag as Tag
from mmalignments.services.io import (
    ensure,
    open_fastq,
    write_fastq_check_results,
    write_json,
)
from mmalignments.models.resources import ResourceConfig  # type: ignore[import]

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
    cfg: PostQCConfig, targets: Element, baits: Element
) -> tuple[Element, Element]:
    """
    Validate and prepare targets and baits for HS metrics and mosdepth.

    Parameters
    ----------
    cfg : PostQCConfig
        Configuration for post-alignment QC.
    targets : Element
        Target intervals for HS metrics and mosdepth.
    baits : Element
        Bait intervals for HS metrics.

    Returns
    -------
    tuple[Element, Element]
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
        if baits is None:
            logger.warning(
                "bait_intervals and/or target_intervals not provided, "
                "using targets for both. This may not be ideal for HS metrics."
            )
            baits = targets

    if cfg.run_mosdepth and targets is None:
        raise ValueError("run_mosdepth=True requires targets to be provided")

    return targets, baits


########################################################################################
# Elements
########################################################################################


def build_input_check(
    sample: NextGenSampleElement,
    n_records: int = 10000,
    *,
    tag: Tag | ElementTag | None = None,
    outdir: Path | str | None = None,
    filename: Path | str | None = None,
) -> Element:
    """Create Element that performs FASTQ sanity checks."""
    default_tag = ElementTag(
        root=sample.tag.root,
        level=sample.tag.level + 1,
        stage=Stage.QC,
        method=Method.CUSTOM,
        state=State.STAT,
        omics=sample.tag.omics,
        ext="json",
    )
    tag = merge_tag(default_tag, tag) if tag is not None else default_tag
    outdir = Path(outdir or sample.fastq1.parent / "qc" / "input_check")

    json_filename = filename or tag.default_output
    out_json = outdir / json_filename
    out_txt = out_json.with_suffix(".txt")

    def _check_fastq() -> None:
        results: dict[str, Any] = {
            "sample_name": sample.name,
            "pairing": sample.pairing.value,
            "checks": {},
            "errors": [],
            "warnings": [],
        }

        try:
            # R1 scan
            r1 = _scan_fastq(sample.fastq_r1, read_label="R1", n_records=n_records)
            results["checks"]["r1_records_checked"] = r1.records_checked
            results["errors"].extend(r1.errors)
            results["checks"]["r1_format_ok"] = len(r1.errors) == 0

            # R2 scan + pairing check (only if paired)
            if sample.pairing == "paired" and sample.fastq_r2:
                r2 = _scan_fastq(sample.fastq_r2, read_label="R2", n_records=n_records)
                results["checks"]["r2_records_checked"] = r2.records_checked
                results["errors"].extend(r2.errors)
                results["checks"]["r2_format_ok"] = len(r2.errors) == 0

                mismatches, mismatch_warnings = _pairing_mismatches(r1.ids, r2.ids)
                results["checks"]["id_mismatches"] = mismatches
                results["checks"]["ids_match"] = mismatches == 0
                results["warnings"].extend(mismatch_warnings)

                if mismatches > 0:
                    results["warnings"].append(
                        f"Found {mismatches} ID mismatches in {min(r1.records_checked, r2.records_checked)} checked records"  # noqa: E501
                    )

            results["status"] = "PASS" if not results["errors"] else "FAIL"

        except Exception as e:
            results["status"] = "ERROR"
            results["errors"].append(f"Exception during check: {e}")
            logger.error(f"Error checking FASTQ files for {sample.name}: {e}")

        write_json(out_json, results)
        write_fastq_check_results(out_txt, sample.name, sample.pairing, results)

        logger.info(f"Input check completed for {sample.name}: {results['status']}")

    input_files = [sample.fastq_r1] + ([sample.fastq_r2] if sample.fastq_r2 else [])
    _check_fastq.command = ["build_input_check"]
    return Element(
        key=f"{tag.default_name}_build_input_check",
        run=_check_fastq,
        tag=tag,
        determinants=[f"n_records={n_records}"],
        inputs=input_files,
        artifacts={"input_check_json": out_json, "input_check_txt": out_txt},
        pres=[],
    )


###############################################################################
# PreAlignment QC workflow
###############################################################################


def pre_alignment_qc(
    sample: NextGenSampleElement,
    *,
    tags: Mapping[str, Tag] | None = None,
    outdir: Path | str | None = None,
    filename: Path | str | None = None,
    parameters: Mapping[str, Params] | None = None,
    preqc_cfg: PreQCConfig | None = None,
    cfgs: Mapping[str, ExternalRunConfig] | None = None,
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
    sample : Sample
        Sample object containing FASTQ file paths and metadata.
    parameters : Mapping[str, Params] | None
        Dictionary of parameter overrides for each QC step, keyed by step
        name.
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
    >>> from mmalignments.models.data import Sample
    >>> from mmalignments.models.externals import ExternalRunConfig
    >>> sample = Sample(name="sample1", pairing="paired",
    ...                 fastq_r1_path="raw_R1.fastq.gz",
    ...                 fastq_r2_path="raw_R2.fastq.gz")
    >>> qc_elements = pre_alignment_qc(
    ...     sample,
    ...     outdir="results/qc",
    ...     preqc_cfg=PreQCConfig(threads=8),
    ... )
    >>> for elem in qc_elements:
    ...     elem.run()
    """
    preqc_cfg = preqc_cfg or PreQCConfig()
    outdir = outdir or preqc_cfg.qc_root or (sample.result_dir / "qc")
    qc_elements = []
    tags = tags or {}
    parameters = parameters or {}
    cfgs = cfgs or {}

    # 1. Input check
    logger.info(f"Creating input check for {sample.name}")
    if preqc_cfg.run_sanity_check:
        logger.info(f"Adding input sanity check for {sample.name}")
        input_check = build_input_check(
            sample=sample,
            outdir=outdir,
            n_records=preqc_cfg.n_records,
        )
        qc_elements.insert(0, input_check)

    # 2. fastp Element
    if preqc_cfg.run_fastp:
        logger.info(f"Creating fastp QC for {sample.name}")
        fastp = FastP()
        fastp_elem = fastp.qc(
            sample=sample,
            tag=tags.get("fastp", None),
            outdir=outdir,
            params=parameters.get("fastp", Params()),
            cfg=cfgs.get("fastp", ExternalRunConfig(threads=preqc_cfg.threads)),
        )
        qc_elements.insert(0, fastp_elem)

    # 3. FastQC on raw FASTQs
    fastqc = FastQC()
    if preqc_cfg.run_fastqc_raw:
        logger.info(f"Creating FastQC for raw reads of {sample.name}")
        fastqc_raw = fastqc.qc(
            sample=sample,
            label="raw",
            tag=tags.get("fastqc_raw", None),
            outdir=outdir,
            params=parameters.get("fastqc_raw", Params()),
            cfg=cfgs.get("fastqc_raw", ExternalRunConfig(threads=preqc_cfg.threads)),
        )
        qc_elements.insert(0, fastqc_raw)

    # 4. FastQC on cleaned FASTQs
    if preqc_cfg.run_fastqc_cleaned:
        logger.info(f"Creating FastQC for cleaned reads of {sample.name}")
        fastqc_cleaned = fastqc.qc(
            sample=fastp_elem,  # Use fastp element as input (has cleaned FASTQs)
            label="cleaned",
            tag=tags.get("fastqc_cleaned", None),
            outdir=outdir,
            params=parameters.get("fastqc_cleaned", Params()),
            cfg=cfgs.get(
                "fastqc_cleaned", ExternalRunConfig(threads=preqc_cfg.threads)
            ),  # noqa: E501
        )
        qc_elements.insert(0, fastqc_cleaned)

    # 5. QC summary Element (pres=[fastp, fastqc_raw, fastqc_cleaned])
    if preqc_cfg.run_summary:
        logger.info(f"Creating QC summary for {sample.name}")
        qc_summary = build_qc_summary(
            sample=sample,
            elements_to_summary=[fastp_elem, fastqc_raw, fastqc_cleaned],
            outdir=outdir,
        )
        qc_elements.insert(0, qc_summary)

    logger.info(
        f"Created pre-alignment QC workflow with {len(qc_elements)} steps for {sample.name}"  # noqa: E501
    )
    return tuple(qc_elements)


########################################################################################
# PostAlignment QC workflow
########################################################################################


def post_mapping_qc_with_multiqc(
    mapped: MappedElement,
    reference: Genome,
    *,
    refdict_element: Element | None = None,
    targets: Element | Path | str | None = None,
    baits: Element | Path | str | None = None,
    tags: Mapping[str, Tag | ElementTag] | None = None,
    outdir: Path | str | None = None,
    filename: Path | str | None = None,
    main_qc_dir: Path | str | None = None,
    parameters: Mapping[str, Params] | None = None,
    externalruncfgs: Mapping[str, ExternalRunConfig | None] | None = None,
    postqc_cfg: PostQCConfig | None = None,
) -> tuple[Element, ...]:
    """Run post-mapping QC and aggregate results with MultiQC.

    This is a convenience wrapper around post_mapping_qc() that also creates
    a MultiQC Element to aggregate all QC reports into a single HTML report.

    Parameters
    ----------
    mapped : MappedElement
        MappedElement containing the BAM file and metadata.
    reference : Genome
        Reference genome object with FASTA file path.
    targets : Element | Path | str | None
        BED file or Element containing target regions.
    baits : Element | Path | str | None
        Interval list for bait regions (Picard format).
    main_qc_dir : Path | str | None
        Base directory for QC outputs. If None, defaults to BAM parent / "qc".
    outdir : Path | str | None
        Output directory for MultiQC report.
    filename : str | None
        Custom name for the MultiQC report HTML file
        (default: "pre_alignment_multiqc.html").
    parameters : Mapping[str, Params] | None
        Dictionary of parameters for QC steps and MultiQC, keyed by step name.
    postqc_cfg : PostQCConfig | None, optional
        Configuration object for post-mapping QC, by default PostQCConfig()
    multiqc_output_dir : Path | str | None
        Output directory for MultiQC report (default: "results/qc/multiqc").
    multiqc_report_name : str | None
        Custom name for the MultiQC report HTML file.
    parameters : Mapping[str, Params] | None
        Dictionary of parametersfor QC steps and MultiQC, keyed by step name.
        For example:
        {
            "gatk": Params(...),
            "multiqc": Params(...),
        }
    externalruncfgs : Mapping[str, ExternalRunConfig | None] | None
        Dictionary of ExternalRunConfig overrides for each QC step and MultiQC, keyed by
        step name. For example:
        {
            "gatk": ExternalRunConfig(...),
            "multiqc": ExternalRunConfig(...),
        }

    Returns
    -------
    tuple[Element, ...]
            Tuple of (multiqc_element, *qc_elements) where qc_elements is the
            list of individual QC Elements and multiqc_element is the MultiQC
            aggregation Element that depends on all QC elements. The first
            returned element in last Element in the QC chain and depends on the
            previous elements.

    Examples
    --------
    >>> qc_elems, multiqc_elem = post_mapping_qc_with_multiqc(
    ...     mapped=mapped_sample,
    ...     reference=genome,
    ...     targets=targets_bed,
    ...     baits=targets_interval_list,
    ...     threads_mosdepth=8
    ... )
    >>> # Execute all QC steps
    >>> for elem in qc_elems:
    ...     elem.run()
    >>> # Generate MultiQC report
    >>> multiqc_elem.run()
    """
    postqc_cfg = postqc_cfg or PostQCConfig(qc_root=main_qc_dir)
    analysis_dir = main_qc_dir or postqc_cfg.qc_root or mapped.bam.parent / "qc"
    parameters = parameters or {}
    tag = tags.pop("multiqc", None) if tags else None
    cfg = (
        externalruncfgs.get("multiqc", ExternalRunConfig())
        if externalruncfgs
        else ExternalRunConfig()
    )
    params = (
        parameters.pop("multiqc", Params(force=True))
        if parameters
        else Params(force=True)
    )
    qc_elements = post_mapping_qc(
        mapped=mapped,
        reference=reference,
        refdict_element=refdict_element,
        targets=targets,
        baits=baits,
        tags=tags,
        outdir=main_qc_dir,
        parameters=parameters,
        externalruncfgs=externalruncfgs,
        postqc_cfg=postqc_cfg,
    )

    # Create MultiQC aggregation element
    multiqc = MultiQC()
    params = parameters.get("multiqc", Params()) if parameters else Params()
    multiqc_elem = multiqc.aggregate(
        analysis_dir=analysis_dir,
        qc_elements=qc_elements,
        tag=tag,
        outdir=outdir,
        filename=filename,
        params=params,
        cfg=cfg,
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
    refdict_element: Element | None = None,
    targets: Element | Path | str | None = None,
    baits: Element | Path | str | None = None,
    tags: Mapping[str, Tag] | None = None,
    outdir: Path | str | None = None,
    parameters: Mapping[str, Params] | None = None,
    externalruncfgs: Mapping[str, ExternalRunConfig | None] | None = None,
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
        Run Picard CollectAlignmentSummaryMetrics (default: True).
    tags : Optional[Mapping[str, Tag]], optional
        Dictionary of Tag overrides for each QC step, keyed by step name.
    outdir : Optional[Path | str], optional
        Base directory for QC outputs. If None, defaults to BAM parent / "qc".
    filename : Optional[Path | str], optional
        Base filename for QC outputs. If None, defaults to using BAM stem with
        appropriate suffixes for each QC step.
    parameters : Optional[Mapping[str, Params]]
        Dictionary of parameter overrides for each QC step, keyed by step name.
        For example:
        {
            "gatk": Params(...),
            "mosdepth": Params(...),
            "samtools": Params(...),
        }
    externalruncfgs : Optional[Mapping[str, ExternalRunConfig | None]]
        Dictionary of ExternalRunConfig overrides for each QC step, keyed by step
        name.
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
    parameters = parameters or {}
    external_cfgs = externalruncfgs or {}
    tags = tags or {}
    output_dir = outdir or postqc_cfg.qc_root or (mapped.bam.parent / "qc")
    picard_dir = output_dir / "picard"
    mosdepth_dir = output_dir / "mosdepth"
    samtools_dir = output_dir / "samtools"
    # threads_mosdepth = cfg.threads_mosdepth or 1
    # no_per_base = cfg.no_per_base or True

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
            tag=tags.get("alignmetrics", None),
            outdir=picard_dir,
            params=parameters.get("alignmetrics", Params()),
            cfg=external_cfgs.get("alignmetrics", ExternalRunConfig()),
        )
        qc_elements["alignmetrics"] = alignment_elem

    # 2. Picard CollectInsertSizeMetrics
    if postqc_cfg.run_insert_size:
        logger.info(f"Creating insert size metrics QC for {mapped.name}")
        insert_elem = gatk.insertmetrics(
            mapped=mapped,
            tag=tags.get("insertmetrics", None),
            outdir=picard_dir,
            params=parameters.get("insertmetrics", Params()),
            cfg=external_cfgs.get("insertmetrics", ExternalRunConfig()),
        )
        qc_elements["insertmetrics"] = insert_elem

    # 3. Picard CollectHsMetrics
    if postqc_cfg.run_hs_metrics:
        # create a sequence dict ... once!
        if refdict_element is None:
            raise ValueError(
                "run_hs_metrics=True requires a sequence dictionary for the reference. "
                "Please provide refdict_element or create one using GATK's CreateSequenceDictionary."
            )
        # create intervals from bed
        targets_intervals = gatk.bed2interval(
            bed_element=targets,
            sequence_dict=refdict_element,
            tag=tags.get("targets2interval", None),
            params=parameters.get("targets2interval", Params()),
            cfg=external_cfgs.get("targets2interval", ExternalRunConfig()),
        )
        if baits is not None:
            # create another interval list for baits (Picard format)
            bait_intervals = gatk.bed2interval(
                bed_element=baits,
                sequence_dict=refdict_element,
                tag=tags.get("targets2interval", None),
                params=parameters.get("targets2interval", Params()),
                cfg=external_cfgs.get("targets2interval", ExternalRunConfig()),
            )
        else:
            bait_intervals = None
        logger.info(f"Creating HS metrics QC for {mapped.name}")
        hs_elem = gatk.hs(
            mapped=mapped,
            reference=reference,
            baits=bait_intervals,
            targets=targets_intervals,
            tag=tags.get("hs", None),
            outdir=picard_dir,
            params=parameters.get("hs", Params()),
            cfg=external_cfgs.get("hs", ExternalRunConfig()),
        )
        qc_elements["hs"] = hs_elem

    # 4. mosdepth coverage
    if postqc_cfg.run_mosdepth:
        logger.info(f"Creating mosdepth coverage QC for {mapped.name}")
        params = parameters.get("mosdepth", Params(threads=47, no_per_base=True))
        if targets and not "by" in params:
            params = params.override(by=targets.bed)
        depth_elem = mosdepth.coverage(
            mapped=mapped,
            targets=targets,
            tag=tags.get("mosdepth", None),
            outdir=mosdepth_dir,  # Use default path
            params=parameters.get("mosdepth", params),
            cfg=external_cfgs.get("mosdepth", ExternalRunConfig()),
        )
        print(depth_elem.run.command_display)
        qc_elements["mosdepth"] = depth_elem

    # 5. samtools flagstat
    if postqc_cfg.run_samtools_flagstat:
        logger.info(f"Creating samtools flagstat QC for {mapped.name}")
        flagstat_elem = st.flagstat(
            mapped=mapped,
            tag=tags.get("flagstat", None),
            outdir=samtools_dir,  # Use default path
            params=parameters.get("flagstat", Params()),
            cfg=external_cfgs.get("flagstat", ExternalRunConfig()),
        )
        qc_elements["flagstat"] = flagstat_elem

    # 6. samtools stats
    if postqc_cfg.run_samtools_stats:
        logger.info(f"Creating samtools stats QC for {mapped.name}")
        stats_elem = st.stats(
            mapped=mapped,
            outdir=samtools_dir,
            tag=tags.get("stats", None),
            params=parameters.get("stats", Params()),
            cfg=external_cfgs.get("stats", ExternalRunConfig()),
        )
        qc_elements["stats"] = stats_elem

    logger.info(f"Created {len(qc_elements)} QC elements for {mapped.name}")
    return qc_elements
