"""Parsers for QC tool outputs.

This module provides functions to parse QC tool outputs (fastp JSON, FastQC ZIP)
and extract relevant metrics in a consistent format for QC summary generation.
"""

from __future__ import annotations

import json
import logging
import zipfile
from inspect import signature
from pathlib import Path
from typing import Any, TextIO

from mmalignments.models.data import Pairing
from mmalignments.models.elements import Element, NextGenSampleElement
from mmalignments.models.tags import (
    ElementTag,
    Method,
    PartialElementTag,
    Stage,
    State,
    merge_tag,
)

logger = logging.getLogger(__name__)


def parse_fastq_record(fh: TextIO) -> tuple[str, str, str, str] | None:
    """Parse one FASTQ record (4 lines).

    Parameters
    ----------
    fh : TextIO
        Open file-handle positioned at the start of a FASTQ record.

    Returns
    -------
    tuple[str, str, str, str] | None
        A tuple of `(header, seq, plus, qual)` for the record, or ``None``
        when end-of-file is reached.
    """
    header = fh.readline().strip()
    if not header:
        return None
    seq = fh.readline().strip()
    plus = fh.readline().strip()
    qual = fh.readline().strip()
    return (header, seq, plus, qual)


def parse_fastp_json(json_path: Path | str) -> dict[str, Any]:
    """Parse fastp JSON output and extract key QC metrics.

    Parameters
    ----------
    json_path : Path | str
        Path to the fastp JSON output file.

    Returns
    -------
    dict[str, Any]
        Dictionary containing extracted metrics, e.g.:
        - total_reads_before, total_reads_after
        - total_bases_before, total_bases_after
        - q20_rate_before/after, q30_rate_before/after
        - gc_content_before/after
        - duplication_rate
        - adapter_trimmed_reads, adapter_trimmed_bases

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


def parse_fastqc_zip(zip_path: Path | str) -> dict[str, Any]:
    """Parse FastQC ZIP output and extract key QC metrics from summary.txt.

    Parameters
    ----------
    zip_path : Path | str
        Path to the FastQC ZIP file (e.g., sample_fastqc.zip).

    Returns
    -------
    dict[str, Any]
        Dictionary containing extracted metrics, for example module PASS/WARN/FAIL
        statuses from ``summary.txt`` and basic statistics such as
        ``total_sequences``, ``sequence_length`` and ``gc_content`` from
        ``fastqc_data.txt`` when available.

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


def write_build_qc_tsv(
    out_tsv: Path | str, sample: NextGenSampleElement, summary: dict[str, Any]
) -> Path:
    """Write a flattened TSV summary from the QC summary dict.

    Parameters
    ----------
    out_tsv : Path | str
        Path to the output TSV file to write.
    sample : NextGenSampleElement
        Sample element providing sample metadata (name, pairing, ...).
    summary : dict[str, Any]
        QC summary dictionary as produced by :func:`build_qc_summary`.

    Returns
    -------
    Path
        The path to the written TSV file.
    """
    out_tsv = Path(out_tsv)
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
    return out_tsv


def build_qc_summary(
    sample: NextGenSampleElement,
    elements_to_summary: tuple[Element, ...],
    *,
    tag: PartialElementTag | ElementTag | None = None,
    outdir: Path | str | None,
    filename: str | None = None,
) -> Element:
    """Build QC summary from fastp and FastQC outputs (pure Python Element).

    Creates an Element that parses QC tool outputs and generates summary
    files in JSON and TSV formats.

    Parameters
    ----------
    sample : NextGenSampleElement
        NextGenSampleElement object representing a sample and containing metadata.
    elements_to_summary : tuple[Element, ...]
        Tuple of QC Elements whose outputs should be included in the QC summary.
    tag : Tag | ElementTag | None
        Partial or full Element tag for the output Element, used for default
        naming. If not provided, a default tag will be generated based on
        the input Element's name and level.
    outdir : Path | str | None
        Directory for the sorted BAM file. If not provided, defaults to
        the same directory as the input BAM file.
    filename : Path | str | None
        Filename override. If not provided, defaults to ``tag.default_output``.

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
    default_tag = ElementTag(
        root=sample.tag.root,
        level=max([el.tag.level for el in elements_to_summary]) + 1,
        stage=Stage.QC,
        method=Method.MULTIQC,
        state=State.REPORT,
        omics=sample.tag.omics,
        ext="json",
    )
    tag = merge_tag(default_tag, tag) if tag is not None else default_tag
    outdir = Path(outdir) or sample.result_dir / "qc"
    out_json = outdir / (filename or (tag.default_name + ".json"))
    out_tsv = out_json.with_suffix(".tsv")
    input_files = []
    for qc_element in elements_to_summary:
        input_files.extend(qc_element.output_files)

    def _build_summary() -> None:
        """Generate QC summary files."""
        summary = {
            "sample_name": sample.name,
            "pairing": sample.pairing.value,
        }
        for qc_element in elements_to_summary:
            # parse FastP json
            method = qc_element.tag.method
            if method == Method.FASTP:
                fastp_metrics = parse_fastp_json(qc_element.json)
                summary["fastp"] = fastp_metrics

            elif method == Method.FASTQC:
                # # Parse FastQC ZIPs
                label = qc_element.tag.param
                summary[f"fastqc_{label}_r1"] = parse_fastqc_zip(
                    qc_element.artifacts["zip_r1"]
                )
                if sample.pairing == Pairing.PAIRED:
                    summary[f"fastqc_{label}_r2"] = parse_fastqc_zip(
                        qc_element.artifacts["zip_r2"]
                    )
            else:
                raise NotImplementedError(f"Unsupported QC method: {method}")
        # Write JSON summary
        with open(out_json, "w") as f:
            json.dump(summary, f, indent=2)

        # Write TSV summary (flattened)
        write_build_qc_tsv(out_tsv, sample, summary)
        logger.info(f"QC summary written to {out_json} and {out_tsv}")

    _build_summary.command = [str(signature(_build_summary))]  # type: ignore[attr-defined]
    _build_summary.command_display = f"Build QC summary for {sample.name}"
    key = f"{tag.default_name}_qc_summary"
    return Element(
        key=key,
        run=_build_summary,
        tag=tag,
        inputs=input_files,
        artifacts={
            "qc_summary_json": out_json,
            "qc_summary_tsv": out_tsv,
        },
        pres=elements_to_summary,
    )
