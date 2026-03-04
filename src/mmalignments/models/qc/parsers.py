"""Parsers for QC tool outputs.

This module provides functions to parse QC tool outputs (fastp JSON, FastQC ZIP)
and extract relevant metrics in a consistent format for QC summary generation.
"""

from __future__ import annotations

import json
import logging
import zipfile
from pathlib import Path
from typing import Any, Dict, Optional, Union
from inspect import signature
from mmalignments.models.data import Sample
from mmalignments.models.tasks import Element
from mmalignments.services.io import ensure

logger = logging.getLogger(__name__)


def parse_fastq_record(fh):
    """Parse one FASTQ record (4 lines)."""
    header = fh.readline().strip()
    if not header:
        return None
    seq = fh.readline().strip()
    plus = fh.readline().strip()
    qual = fh.readline().strip()
    return (header, seq, plus, qual)


def parse_fastp_json(json_path: Union[Path, str]) -> Dict[str, Any]:
    """Parse fastp JSON output and extract key QC metrics.

    Parameters
    ----------
    json_path : Union[Path, str]
            Path to fastp JSON output file.

    Returns
    -------
    Dict[str, Any]
            Dictionary containing extracted metrics:
            - total_reads_before: Total reads before filtering
            - total_reads_after: Total reads after filtering
            - total_bases_before: Total bases before filtering
            - total_bases_after: Total bases after filtering
            - q20_rate_before: Q20 rate before filtering
            - q20_rate_after: Q20 rate after filtering
            - q30_rate_before: Q30 rate before filtering
            - q30_rate_after: Q30 rate after filtering
            - gc_content_before: GC content before filtering
            - gc_content_after: GC content after filtering
            - duplication_rate: Duplication rate
            - adapter_trimmed_reads: Number of reads with adapters trimmed
            - adapter_trimmed_bases: Number of bases trimmed from adapters

    Examples
    --------
    >>> metrics = parse_fastp_json("sample_fastp.json")
    >>> print(f"Total reads after: {metrics['total_reads_after']}")
    >>> print(f"Q30 rate after: {metrics['q30_rate_after']:.2%}")
    """
    json_path = Path(json_path)

    if not json_path.exists():
        logger.warning(f"fastp JSON file not found: {json_path}")
        return {}

    try:
        with open(json_path, "r") as f:
            data = json.load(f)

        metrics = {}

        # Summary metrics
        if "summary" in data:
            summary = data["summary"]
            if "before_filtering" in summary:
                bf = summary["before_filtering"]
                metrics["total_reads_before"] = bf.get("total_reads", 0)
                metrics["total_bases_before"] = bf.get("total_bases", 0)
                metrics["q20_rate_before"] = bf.get("q20_rate", 0)
                metrics["q30_rate_before"] = bf.get("q30_rate", 0)
                metrics["gc_content_before"] = bf.get("gc_content", 0)

            if "after_filtering" in summary:
                af = summary["after_filtering"]
                metrics["total_reads_after"] = af.get("total_reads", 0)
                metrics["total_bases_after"] = af.get("total_bases", 0)
                metrics["q20_rate_after"] = af.get("q20_rate", 0)
                metrics["q30_rate_after"] = af.get("q30_rate", 0)
                metrics["gc_content_after"] = af.get("gc_content", 0)

        # Duplication metrics
        if "duplication" in data:
            metrics["duplication_rate"] = data["duplication"].get("rate", 0)

        # Adapter trimming metrics
        if "adapter_cutting" in data:
            ac = data["adapter_cutting"]
            metrics["adapter_trimmed_reads"] = ac.get("adapter_trimmed_reads", 0)
            metrics["adapter_trimmed_bases"] = ac.get("adapter_trimmed_bases", 0)

        # Filtering result
        if "filtering_result" in data:
            fr = data["filtering_result"]
            metrics["passed_filter_reads"] = fr.get("passed_filter_reads", 0)
            metrics["low_quality_reads"] = fr.get("low_quality_reads", 0)
            metrics["too_many_N_reads"] = fr.get("too_many_N_reads", 0)
            metrics["too_short_reads"] = fr.get("too_short_reads", 0)

        return metrics

    except Exception as e:
        logger.error(f"Error parsing fastp JSON {json_path}: {e}")
        return {}


def parse_fastqc_zip(zip_path: Union[Path, str]) -> Dict[str, Any]:
    """Parse FastQC ZIP output and extract key QC metrics from summary.txt.

    Parameters
    ----------
    zip_path : Union[Path, str]
            Path to FastQC ZIP file (e.g., sample_fastqc.zip).

    Returns
    -------
    Dict[str, Any]
            Dictionary containing extracted metrics:
            - basic_statistics: PASS/WARN/FAIL status
            - per_base_sequence_quality: PASS/WARN/FAIL status
            - per_sequence_quality_scores: PASS/WARN/FAIL status
            - per_base_sequence_content: PASS/WARN/FAIL status
            - per_sequence_gc_content: PASS/WARN/FAIL status
            - per_base_n_content: PASS/WARN/FAIL status
            - sequence_length_distribution: PASS/WARN/FAIL status
            - sequence_duplication_levels: PASS/WARN/FAIL status
            - overrepresented_sequences: PASS/WARN/FAIL status
            - adapter_content: PASS/WARN/FAIL status
            - total_sequences: Total number of sequences
            - sequence_length: Sequence length range
            - gc_content: %GC content

    Examples
    --------
    >>> metrics = parse_fastqc_zip("sample_R1_fastqc.zip")
    >>> print(f"Per base quality: {metrics['per_base_sequence_quality']}")
    >>> print(f"Total sequences: {metrics['total_sequences']}")
    """
    zip_path = Path(zip_path)

    if not zip_path.exists():
        logger.warning(f"FastQC ZIP file not found: {zip_path}")
        return {}

    try:
        metrics = {}

        with zipfile.ZipFile(zip_path, "r") as zf:
            # Find summary.txt in the ZIP
            summary_file = None
            for name in zf.namelist():
                if name.endswith("summary.txt"):
                    summary_file = name
                    break

            if summary_file is None:
                logger.warning(f"summary.txt not found in {zip_path}")
                return {}

            # Parse summary.txt
            with zf.open(summary_file) as f:
                for line in f:
                    line = line.decode("utf-8").strip()
                    if not line:
                        continue

                    parts = line.split("\t")
                    if len(parts) >= 2:
                        status = parts[0]  # PASS, WARN, or FAIL
                        module_name = parts[1]

                        # Convert module name to metric key
                        key = module_name.lower().replace(" ", "_")
                        metrics[key] = status

            # Parse fastqc_data.txt for additional metrics
            data_file = None
            for name in zf.namelist():
                if name.endswith("fastqc_data.txt"):
                    data_file = name
                    break

            if data_file:
                with zf.open(data_file) as f:
                    in_basic_stats = False
                    for line in f:
                        line = line.decode("utf-8").strip()

                        if line.startswith(">>Basic Statistics"):
                            in_basic_stats = True
                            continue
                        elif line.startswith(">>END_MODULE"):
                            in_basic_stats = False
                            continue

                        if in_basic_stats and line and not line.startswith("#"):
                            parts = line.split("\t")
                            if len(parts) >= 2:
                                measure = parts[0]
                                value = parts[1]

                                if measure == "Total Sequences":
                                    metrics["total_sequences"] = int(value)
                                elif measure == "Sequence length":
                                    metrics["sequence_length"] = value
                                elif measure == "%GC":
                                    metrics["gc_content"] = float(value)

        return metrics

    except Exception as e:
        logger.error(f"Error parsing FastQC ZIP {zip_path}: {e}")
        return {}


def build_qc_summary(
    sample: Sample,
    fastp_json: Optional[Union[Path, str]],
    fastqc_raw_r1_zip: Optional[Union[Path, str]],
    fastqc_raw_r2_zip: Optional[Union[Path, str]],
    fastqc_cleaned_r1_zip: Optional[Union[Path, str]],
    fastqc_cleaned_r2_zip: Optional[Union[Path, str]],
    out_dir: Union[Path, str],
    out_json: str,
    out_tsv: str,
) -> Element:
    """Build QC summary from fastp and FastQC outputs (pure Python Element).

    Creates an Element that parses QC tool outputs and generates summary
    files in JSON and TSV formats.

    Parameters
    ----------
    sample : Sample
            Sample object containing metadata.
    fastp_json : Optional[Union[Path, str]]
            Path to fastp JSON output.
    fastqc_raw_r1_zip : Optional[Union[Path, str]]
            Path to FastQC ZIP for raw R1.
    fastqc_raw_r2_zip : Optional[Union[Path, str]]
            Path to FastQC ZIP for raw R2 (if paired-end).
    fastqc_cleaned_r1_zip : Optional[Union[Path, str]]
            Path to FastQC ZIP for cleaned R1.
    fastqc_cleaned_r2_zip : Optional[Union[Path, str]]
            Path to FastQC ZIP for cleaned R2 (if paired-end).
    out_dir : Union[Path, str]
            Output directory for summary files.
    out_json : Union[Path, str]
            Output JSON summary file path.
    out_tsv : Union[Path, str]
            Output TSV summary file path.

    Returns
    -------
    Element
            An Element that generates QC summary when run. Artifacts include:
            - qc_summary_json: JSON summary file
            - qc_summary_tsv: TSV summary file

    Examples
    --------
    >>> summary_elem = build_qc_summary(
    ...     sample=sample,
    ...     fastp_json="fastp.json",
    ...     fastqc_raw_r1_zip="raw_R1_fastqc.zip",
    ...     fastqc_raw_r2_zip="raw_R2_fastqc.zip",
    ...     fastqc_cleaned_r1_zip="cleaned_R1_fastqc.zip",
    ...     fastqc_cleaned_r2_zip="cleaned_R2_fastqc.zip",
    ...     out_json="sample_qc.json",
    ...     out_tsv="sample_qc.tsv"
    ... )
    >>> summary_elem.run()
    """
    ensure(out_dir)
    out_json = (Path(out_dir) / out_json).absolute()
    out_tsv = (Path(out_dir) / out_tsv).absolute()

    # Collect input files
    input_files = []
    if fastp_json:
        input_files.append(Path(fastp_json))
    if fastqc_raw_r1_zip:
        input_files.append(Path(fastqc_raw_r1_zip))
    if fastqc_raw_r2_zip:
        input_files.append(Path(fastqc_raw_r2_zip))
    if fastqc_cleaned_r1_zip:
        input_files.append(Path(fastqc_cleaned_r1_zip))
    if fastqc_cleaned_r2_zip:
        input_files.append(Path(fastqc_cleaned_r2_zip))

    def _build_summary() -> None:
        """Generate QC summary files."""
        summary = {
            "sample_name": sample.name,
            "pairing": sample.pairing,
        }

        # Parse fastp JSON
        if fastp_json:
            fastp_metrics = parse_fastp_json(fastp_json)
            summary["fastp"] = fastp_metrics

        # Parse FastQC ZIPs
        if fastqc_raw_r1_zip:
            summary["fastqc_raw_r1"] = parse_fastqc_zip(fastqc_raw_r1_zip)
        if fastqc_raw_r2_zip:
            summary["fastqc_raw_r2"] = parse_fastqc_zip(fastqc_raw_r2_zip)
        if fastqc_cleaned_r1_zip:
            summary["fastqc_cleaned_r1"] = parse_fastqc_zip(fastqc_cleaned_r1_zip)
        if fastqc_cleaned_r2_zip:
            summary["fastqc_cleaned_r2"] = parse_fastqc_zip(fastqc_cleaned_r2_zip)

        # Write JSON summary
        with open(out_json, "w") as f:
            json.dump(summary, f, indent=2)

        # Write TSV summary (flattened)
        with open(out_tsv, "w") as f:
            f.write("sample\tmetric\tvalue\n")
            f.write(f"{sample.name}\tsample_name\t{sample.name}\n")
            f.write(f"{sample.name}\tpairing\t{sample.pairing}\n")

            # Write fastp metrics
            if "fastp" in summary:
                for key, value in summary["fastp"].items():
                    f.write(f"{sample.name}\tfastp_{key}\t{value}\n")

            # Write FastQC metrics
            for qc_type in [
                "fastqc_raw_r1",
                "fastqc_raw_r2",
                "fastqc_cleaned_r1",
                "fastqc_cleaned_r2",
            ]:
                if qc_type in summary:
                    for key, value in summary[qc_type].items():
                        f.write(f"{sample.name}\t{qc_type}_{key}\t{value}\n")

        logger.info(f"QC summary written to {out_json} and {out_tsv}")

    _build_summary.command = [str(signature(_build_summary))]  # type: ignore[attr-defined]
    _build_summary.command_display = f"Build QC summary for {sample.name}"
    return Element(
        name=f"{sample.name}.QC(summary)",
        key=f"{sample.name}_qc_summary",
        run=_build_summary,
        inputs=input_files,
        artifacts={
            "qc_summary_json": out_json,
            "qc_summary_tsv": out_tsv,
        },
        pres=[],
    )
