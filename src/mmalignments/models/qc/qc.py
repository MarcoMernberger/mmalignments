"""Convenience functions for post-mapping quality control workflows.

This module provides high-level functions that chain multiple QC steps
into a single workflow, similar to the alignsort pattern in BWAMem2.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from mmalignments.models.aligners.samtools import Samtools
from mmalignments.models.callers.gatk import GATK
from mmalignments.models.data import Genome, Sample
from mmalignments.models.externals import ExternalRunConfig
from mmalignments.models.parameters import Params
from mmalignments.models.qc.fastp import FastP
from mmalignments.models.qc.fastqc import FastQC
from mmalignments.models.qc.mosdepth import Mosdepth
from mmalignments.models.qc.multiqc import MultiQC
from mmalignments.models.qc.parsers import build_qc_summary
from mmalignments.models.tasks import Element, MappedElement
from mmalignments.services.io import (
    ensure,
    open_fastq,
    write_fastq_check_results,
    write_json,
)

from .parsers import parse_fastq_record

logger = logging.getLogger(__name__)


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
    """Validate a single FASTQ record and return error messages (if any)."""
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
    """Scan up to n_records from a FASTQ and return IDs, record count, and errors."""
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
    """Compare paired IDs and return mismatch count and warning examples."""
    mismatches = 0
    warnings: list[str] = []

    for i, (id1, id2) in enumerate(zip(r1_ids, r2_ids)):
        if id1 != id2:
            mismatches += 1
            if mismatches <= max_examples:
                warnings.append(f"ID mismatch at record {i+1}: R1='{id1}' R2='{id2}'")

    return mismatches, warnings


def build_input_check(
    sample: "Sample",
    out_dir: Path | str | None = None,
    n_records: int = 10000,
) -> "Element":
    """Create Element that performs FASTQ sanity checks."""
    out_dir_path = (
        (sample.result_dir / "qc" / "input_check") if out_dir is None else Path(out_dir)
    )
    out_dir_path = out_dir_path.absolute()
    ensure(out_dir_path)

    out_json = out_dir_path / f"{sample.name}_input_check.json"
    out_txt = out_dir_path / f"{sample.name}_input_check.txt"

    def _check_fastq() -> None:
        results: dict[str, Any] = {
            "sample_name": sample.name,
            "pairing": sample.pairing,
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
        name=f"{sample.name}.QC(sanity)",
        key=f"{sample.name}_build_input_check",
        run=_check_fastq,
        determinants=[f"n_records={n_records}"],
        inputs=input_files,
        artifacts={"input_check_json": out_json, "input_check_txt": out_txt},
        pres=[],
        store_attributes={"multi": False},
    )


###############################################################################
# PreAlignment QC workflow
###############################################################################


@dataclass(frozen=True)
class PreQCConfig:
    qc_root: Path | None = None
    threads: int = 47
    n_records: int = 10_000

    run_sanity_check: bool = True
    run_fastp: bool = True
    run_fastqc_raw: bool = True
    run_fastqc_cleaned: bool = True
    run_summary: bool = True


def pre_alignment_qc(
    sample: Sample,
    parameters: Mapping[str, Params] | None = None,
    cfg: PreQCConfig | None = None,
) -> tuple[Element]:
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
    parameters : Mapping[str, Params]
        Dictionary of parameter overrides for each QC step, keyed by step
        name.
    cfg : PreQCConfig | None
        Configuration object for pre-alignment QC. If None, defaults will
        be used.

    Returns
    -------
    tuple[Element]
        tuple of QC Elements in dependency order:
        input_check, fastp, fastqc_raw, fastqc_cleaned, qc_summary
        the first element is the last in the dependency chain (qc_summary
        depends on all previous steps).


    Examples
    --------
    >>> from mmalignments.models.data import Sample
    >>> sample = Sample(name="sample1", pairing="paired",
    ...                 fastq_r1_path="raw_R1.fastq.gz",
    ...                 fastq_r2_path="raw_R2.fastq.gz")
    >>> qc_elements = pre_alignment_qc(
    ...     sample, "results/qc", PreQCConfig(threads=8)
    ... )
    >>> for elem in qc_elements:
    ...     elem.run()
    """
    cfg = cfg or PreQCConfig()
    if cfg.qc_root:
        qc_root = Path(cfg.qc_root).absolute()
        sample_qc_dir = qc_root / sample.name
        ensure(sample_qc_dir)
    else:
        sample_qc_dir = sample.result_dir / "qc"

    qc_elements = []

    # 1. Input check
    logger.info(f"Creating input check for {sample.name}")
    if cfg.run_sanity_check:
        logger.info(f"Adding input sanity check for {sample.name}")
        input_check = build_input_check(
            sample=sample,
            out_dir=sample_qc_dir,
            n_records=cfg.n_records,
        )
        qc_elements.insert(0, input_check)
    # 2. fastp Element
    if cfg.run_fastp:
        logger.info(f"Creating fastp QC for {sample.name}")
        fastp = FastP()
        params = parameters.get("fastp", Params()) if parameters else Params()
        fastp_elem = fastp.qc(
            sample=sample,
            folder=sample_qc_dir,
            params=params,
            cfg=ExternalRunConfig(threads=cfg.threads),
        )
        qc_elements.insert(0, fastp_elem)

    # 3. FastQC on raw FASTQs
    fastqc = FastQC()
    params = parameters.get("fastqc", Params()) if parameters else Params()
    if cfg.run_fastqc_raw:
        logger.info(f"Creating FastQC for raw reads of {sample.name}")
        fastqc_raw = fastqc.qc(
            sample=sample,
            folder=sample_qc_dir,
            label="raw",
            params=params,
            cfg=ExternalRunConfig(threads=cfg.threads),
        )
        qc_elements.insert(0, fastqc_raw)

    # 4. FastQC on cleaned FASTQs
    if cfg.run_fastqc_cleaned:
        logger.info(f"Creating FastQC for cleaned reads of {sample.name}")
        fastqc_cleaned = fastqc.qc(
            sample=fastp_elem,  # Use fastp element as input (has cleaned FASTQs)
            folder=sample_qc_dir,
            label="cleaned",
            params=params,
            cfg=ExternalRunConfig(threads=cfg.threads),
        )
        qc_elements.insert(0, fastqc_cleaned)

    # 5. QC summary Element (pres=[fastp, fastqc_raw, fastqc_cleaned])
    if cfg.run_summary:
        logger.info(f"Creating QC summary for {sample.name}")
        params = parameters.get("qc_summary", Params()) if parameters else Params()
        qc_summary = build_qc_summary(
            sample=sample,
            fastp_json=fastp_elem.json,
            fastqc_raw_r1_zip=fastqc_raw.zip_r1,
            fastqc_raw_r2_zip=fastqc_raw.zip_r2,
            fastqc_cleaned_r1_zip=fastqc_cleaned.zip_r1,
            fastqc_cleaned_r2_zip=fastqc_cleaned.zip_r2,
            out_dir=sample_qc_dir or sample.result_dir / "qc",
            out_json=f"{sample.name}_qc.json",
            out_tsv=f"{sample.name}_qc.tsv",
        )
        qc_elements.insert(0, qc_summary)

    logger.info(
        f"Created pre-alignment QC workflow with {len(qc_elements)} steps for {sample.name}"  # noqa: E501
    )
    return tuple(qc_elements)


@dataclass(frozen=True)
class PostQCConfig:
    qc_root: Path | None = None
    run_alignment_summary: bool = True
    run_insert_size: bool = True
    run_hs_metrics: bool = True
    run_mosdepth: bool = True
    run_samtools_flagstat: bool = True
    run_samtools_stats: bool = True
    multiqc: bool = True


def post_mapping_qc_with_multiqc(
    mapped: MappedElement,
    reference: Genome,
    *,
    targets: Element | Path | str | None = None,
    baits: Element | Path | str | None = None,
    main_qc_dir: Path | str | None = None,
    multiqc_output_dir: Path | str | None = None,
    multiqc_report_name: str | None = "pre_alignment_multiqc.html",
    parameters: Mapping[str, Params] | None = None,
    externalruncfgs: Mapping[str, ExternalRunConfig | None] | None = None,
    cfg: PostQCConfig | None = None,
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
    multiqc_output_dir : Path | str | None
        Output directory for MultiQC report (default: "results/qc/multiqc").
    multiqc_report_name : str | None
        Custom name for the MultiQC report HTML file
        (default: "pre_alignment_multiqc.html").
    parameters : Mapping[str, Params] | None
        Dictionary of parameters for QC steps and MultiQC, keyed by step name.
    cfg : PostQCConfig | None, optional
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
    main_qc_dir = mapped.bam.parent / "qc"
    cfg = cfg or PostQCConfig(qc_root=main_qc_dir)
    parameters = parameters or {}
    qc_elements = post_mapping_qc(
        mapped=mapped,
        reference=reference,
        targets=targets,
        baits=baits,
        parameters=parameters,
        externalruncfgs=externalruncfgs,
        cfg=cfg,
    )
    # Create MultiQC aggregation element
    multiqc = MultiQC()
    cfg = (
        externalruncfgs.get("multiqc", ExternalRunConfig())
        if externalruncfgs
        else ExternalRunConfig()
    )
    cfg.qc_root = main_qc_dir
    params = parameters.get("multiqc", Params()) if parameters else Params()
    analysis_dir = main_qc_dir
    multiqc_elem = multiqc.aggregate(
        mapped_element_or_analysis_dir=analysis_dir,
        qc_elements=qc_elements,
        output_dir=multiqc_output_dir or Path("results/qc/multiqc"),
        report_name=multiqc_report_name,
        force=True,
        params=params,
        cfg=cfg,
    )

    logger.info(
        f"Created QC workflow with {len(qc_elements)} QC steps "
        f"and MultiQC aggregation for {mapped.name}"
    )

    return multiqc_elem, *qc_elements


# post alignment qc
def post_mapping_qc(
    mapped: MappedElement,
    reference: Genome,
    *,
    targets: Element | Path | str | None = None,
    baits: Element | Path | str | None = None,
    parameters: Mapping[str, Params] | None = None,
    externalruncfgs: Mapping[str, ExternalRunConfig | None] | None = None,
    cfg: PostQCConfig | None = None,
) -> tuple[Element, ...]:
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
    parameters : Optional[Mapping[str, Params]]
            Dictionary of parameter overrides for each QC step, keyed by step name.
            For example:
            {
                "gatk": Params(...),
                "mosdepth": Params(...),
                "samtools": Params(...),
            }
    exterbnalruncfgs : Optional[Mapping[str, ExternalRunConfig | None]]
            Dictionary of ExternalRunConfig overrides for each QC step, keyed by step name.

    cfg : Optional[PostQCConfig], optional
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
    cfg = cfg or PostQCConfig()
    parameters = parameters or {}
    external_cfgs = externalruncfgs or {}
    output_dir = cfg.qc_root or (mapped.bam.parent / "qc")
    picard_dir = output_dir / "picard"
    mosdepth_dir = output_dir / "mosdepth"
    samtools_dir = output_dir / "samtools"
    ensure(picard_dir, mosdepth_dir, samtools_dir)
    # threads_mosdepth = cfg.threads_mosdepth or 1
    # no_per_base = cfg.no_per_base or True

    qc_elements = []

    # tools
    gatk = GATK()
    mosdepth = Mosdepth()
    st = Samtools()

    # Validate baits/targets for HS metrics and mosdepth
    targets, baits = _validate_targets_baits(cfg, targets, baits)

    # 1. Picard CollectAlignmentSummaryMetrics
    if cfg.run_alignment_summary:
        logger.info(f"Creating alignment summary metrics QC for {mapped.name}")
        alignment_elem = gatk.alignment_summary(
            mapped=mapped,
            reference=reference,
            output_metrics=(
                picard_dir / f"{mapped.bam.stem}_alignment_summary.txt"
            ).absolute(),
            params=parameters.get("alignment_summary", Params()),
            cfg=external_cfgs.get("alignment_summary"),
        )
        qc_elements.insert(0, alignment_elem)

    # 2. Picard CollectInsertSizeMetrics
    if cfg.run_insert_size:
        logger.info(f"Creating insert size metrics QC for {mapped.name}")
        insert_elem = gatk.insertsize(
            mapped=mapped,
            output_metrics=(
                picard_dir / f"{mapped.bam.stem}_insert_size.txt"
            ).absolute(),
            output_histogram=(picard_dir / f"{mapped.bam.stem}_insert_size.pdf"),
            params=parameters.get("insertsize", Params()),
            cfg=external_cfgs.get("insertsize", ExternalRunConfig()),
        )
        qc_elements.insert(0, insert_elem)

    # 3. Picard CollectHsMetrics
    if cfg.run_hs_metrics:
        sequence_dict = gatk.sequence_dict(
            reference=reference,
            output_dict=f"cache/picard/{reference.name}.dict",
        )
        targets_intervals = gatk.bed2interval(
            element=targets,
            sequence_dict=sequence_dict,
            params=parameters.get("bed2interval", Params()),
            cfg=external_cfgs.get("bed2interval", ExternalRunConfig()),
        )
        if baits is not None:
            bait_intervals = gatk.bed2interval(
                element=baits,
                sequence_dict=sequence_dict,
                params=parameters.get("bed2interval", Params()),
                cfg=external_cfgs.get("bed2interval", ExternalRunConfig()),
            )
        else:
            bait_intervals = None
        logger.info(f"Creating HS metrics QC for {mapped.name}")
        hs_elem = gatk.hs(
            mapped=mapped,
            reference=reference,
            baits=bait_intervals,
            targets=targets_intervals,
            output_metrics=(
                picard_dir / f"{mapped.bam.stem}_hs_metrics.txt"
            ).absolute(),  # Use default path
            params=parameters.get("hs_metrics", Params()),
            cfg=external_cfgs.get("hs_metrics", ExternalRunConfig()),
        )
        qc_elements.insert(0, hs_elem)
    # 4. mosdepth coverage
    if cfg.run_mosdepth:
        logger.info(f"Creating mosdepth coverage QC for {mapped.name}")
        depth_elem = mosdepth.coverage(
            mapped=mapped,
            targets=targets,
            output_file_prefix=mosdepth_dir
            / f"{mapped.bam.stem}.targets",  # Use default path
            params=parameters.get("mosdepth", Params(threads=47, no_per_base=True)),
            cfg=external_cfgs.get("mosdepth", ExternalRunConfig()),
        )
        qc_elements.insert(0, depth_elem)

    # 5. samtools flagstat
    if cfg.run_samtools_flagstat:
        logger.info(f"Creating samtools flagstat QC for {mapped.name}")
        flagstat_elem = st.flagstat(
            mapped=mapped,
            output_file=(
                samtools_dir / f"{mapped.bam.stem}_flagstat.txt"
            ).absolute(),  # Use default path
            params=parameters.get("samtools_flagstat", Params()),
            cfg=external_cfgs.get("samtools_flagstat", ExternalRunConfig()),
        )
        qc_elements.insert(0, flagstat_elem)

    # 6. samtools stats
    if cfg.run_samtools_stats:
        logger.info(f"Creating samtools stats QC for {mapped.name}")
        stats_elem = st.stats(
            mapped=mapped,
            output_file=(
                samtools_dir / f"{mapped.bam.stem}_stats.txt"
            ).absolute(),  # Use default path
            params=parameters.get("samtools_stats", Params()),
            cfg=external_cfgs.get("samtools_stats", ExternalRunConfig()),
        )
        qc_elements.insert(0, stats_elem)

    logger.info(f"Created {len(qc_elements)} QC elements for {mapped.name}")
    return qc_elements


def _validate_targets_baits(cfg, targets, baits):
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
