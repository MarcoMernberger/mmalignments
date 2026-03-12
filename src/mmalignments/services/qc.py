"""Quality control services using FastQC and MultiQC."""

import pypipegraph2 as ppg  # type: ignore
from pathlib import Path
from typing import List, Optional, Dict
import subprocess
import json


class QCService:
    """Service for running quality control tools."""

    def __init__(self, threads: int = 4):
        """
        Initialize QC service.

        Args:
            threads: Number of threads for FastQC
        """
        self.threads = threads

    def run_fastqc(
        self,
        fastq_files: List[Path],
        output_dir: Path,
    ) -> ppg.MultiFileGeneratingJob:
        """
        Run FastQC on FASTQ files.

        Args:
            fastq_files: List of FASTQ files to analyze
            output_dir: Output directory for FastQC reports

        Returns:
            pypipegraph job creating HTML and zip files
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # FastQC creates _fastqc.html and _fastqc.zip for each input
        outputs = []
        for fq in fastq_files:
            fq = Path(fq)
            base_name = fq.stem
            if base_name.endswith(".fastq"):
                base_name = base_name[:-6]
            elif base_name.endswith(".fq"):
                base_name = base_name[:-3]

            outputs.append(output_dir / f"{base_name}_fastqc.html")
            outputs.append(output_dir / f"{base_name}_fastqc.zip")

        def run_fastqc_cmd():
            cmd = [
                "fastqc",
                "-t",
                str(self.threads),
                "-o",
                str(output_dir),
            ]
            cmd.extend([str(f) for f in fastq_files])

            subprocess.run(cmd, check=True)

        job = ppg.MultiFileGeneratingJob(outputs, run_fastqc_cmd)
        return job

    def run_fastqc_on_sample(
        self,
        sample_name: str,
        fastq_r1: Optional[Path],
        fastq_r2: Optional[Path],
        output_dir: Path,
    ) -> ppg.MultiFileGeneratingJob:
        """
        Run FastQC on a sample's FASTQ files.

        Args:
            sample_name: Sample name
            fastq_r1: R1 FASTQ file
            fastq_r2: R2 FASTQ file (optional)
            output_dir: Output directory

        Returns:
            pypipegraph job
        """
        fastq_files = []
        if fastq_r1:
            fastq_files.append(fastq_r1)
        if fastq_r2:
            fastq_files.append(fastq_r2)

        sample_dir = output_dir / sample_name
        return self.run_fastqc(fastq_files, sample_dir)

    def run_multiqc(
        self,
        input_dirs: List[Path],
        output_dir: Path,
        report_name: str = "multiqc_report",
        module: Optional[str] = None,
    ) -> ppg.MultiFileGeneratingJob:
        """
        Run MultiQC to aggregate QC reports.

        Args:
            input_dirs: Directories containing QC data (FastQC, Picard, etc.)
            output_dir: Output directory for MultiQC report
            report_name: Name of the report
            module: Specific module to run (e.g., 'fastqc', 'picard')

        Returns:
            pypipegraph job creating MultiQC report
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        report_html = output_dir / f"{report_name}.html"
        report_data = output_dir / f"{report_name}_data"

        outputs = [report_html, report_data]

        def run_multiqc_cmd():
            cmd = [
                "multiqc",
                "-o",
                str(output_dir),
                "-n",
                report_name,
                "--force",
            ]

            if module:
                cmd.extend(["-m", module])

            # Add all input directories
            cmd.extend([str(d) for d in input_dirs])

            subprocess.run(cmd, check=True)

        job = ppg.MultiFileGeneratingJob(outputs, run_multiqc_cmd)
        return job

    def collect_alignment_metrics(
        self,
        bam_file: Path,
        reference_fasta: Path,
        output_file: Path,
    ) -> ppg.FileGeneratingJob:
        """
        Collect alignment summary metrics using Picard.

        Args:
            bam_file: Input BAM file
            reference_fasta: Reference genome
            output_file: Output metrics file

        Returns:
            pypipegraph job
        """
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        def collect_metrics():
            cmd = [
                "picard",
                "CollectAlignmentSummaryMetrics",
                f"INPUT={bam_file}",
                f"OUTPUT={output_file}",
                f"REFERENCE_SEQUENCE={reference_fasta}",
            ]
            subprocess.run(cmd, check=True)

        job = ppg.FileGeneratingJob(output_file, collect_metrics)
        return job

    def collect_insert_size_metrics(
        self,
        bam_file: Path,
        output_file: Path,
        histogram_file: Path,
    ) -> ppg.MultiFileGeneratingJob:
        """
        Collect insert size metrics using Picard.

        Args:
            bam_file: Input BAM file
            output_file: Output metrics file
            histogram_file: Output histogram PDF

        Returns:
            pypipegraph job
        """
        output_file = Path(output_file)
        histogram_file = Path(histogram_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        def collect_metrics():
            cmd = [
                "picard",
                "CollectInsertSizeMetrics",
                f"INPUT={bam_file}",
                f"OUTPUT={output_file}",
                f"HISTOGRAM_FILE={histogram_file}",
            ]
            subprocess.run(cmd, check=True)

        job = ppg.MultiFileGeneratingJob(
            [output_file, histogram_file],
            collect_metrics,
        )
        return job

    def collect_hs_metrics(
        self,
        bam_file: Path,
        reference_fasta: Path,
        bait_intervals: Path,
        target_intervals: Path,
        output_file: Path,
    ) -> ppg.FileGeneratingJob:
        """
        Collect hybrid selection (HS) metrics for exome/targeted sequencing.

        Args:
            bam_file: Input BAM file
            reference_fasta: Reference genome
            bait_intervals: Bait interval list
            target_intervals: Target interval list
            output_file: Output metrics file

        Returns:
            pypipegraph job
        """
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        def collect_metrics():
            cmd = [
                "picard",
                "CollectHsMetrics",
                f"INPUT={bam_file}",
                f"OUTPUT={output_file}",
                f"REFERENCE_SEQUENCE={reference_fasta}",
                f"BAIT_INTERVALS={bait_intervals}",
                f"TARGET_INTERVALS={target_intervals}",
            ]
            subprocess.run(cmd, check=True)

        job = ppg.FileGeneratingJob(output_file, collect_metrics)
        return job

    def run_mosdepth(
        self,
        bam_file: Path,
        output_prefix: Path,
        bed_file: Optional[Path] = None,
        threads: Optional[int] = None,
    ) -> ppg.MultiFileGeneratingJob:
        """
        Run mosdepth for coverage analysis.

        Args:
            bam_file: Input BAM file
            output_prefix: Output prefix (will create multiple files)
            bed_file: BED file for targeted analysis
            threads: Number of threads (default: self.threads)

        Returns:
            pypipegraph job
        """
        if threads is None:
            threads = self.threads

        output_prefix = Path(output_prefix)
        output_prefix.parent.mkdir(parents=True, exist_ok=True)

        # mosdepth creates multiple output files
        outputs = [
            Path(str(output_prefix) + ".mosdepth.global.dist.txt"),
            Path(str(output_prefix) + ".mosdepth.summary.txt"),
        ]

        if bed_file:
            outputs.append(Path(str(output_prefix) + ".regions.bed.gz"))
            outputs.append(Path(str(output_prefix) + ".thresholds.bed.gz"))

        def run_mosdepth_cmd():
            cmd = [
                "mosdepth",
                "-t",
                str(threads),
            ]

            if bed_file:
                cmd.extend(["-b", str(bed_file)])

            cmd.extend(
                [
                    str(output_prefix),
                    str(bam_file),
                ]
            )

            subprocess.run(cmd, check=True)

        job = ppg.MultiFileGeneratingJob(outputs, run_mosdepth_cmd)
        return job

    def comprehensive_bam_qc(
        self,
        sample_name: str,
        bam_file: Path,
        reference_fasta: Path,
        output_dir: Path,
        target_bed: Optional[Path] = None,
        bait_intervals: Optional[Path] = None,
        target_intervals: Optional[Path] = None,
    ) -> Dict[str, ppg.Job]:
        """
        Run comprehensive QC on a BAM file.

        Args:
            sample_name: Sample name
            bam_file: Input BAM file
            reference_fasta: Reference genome
            output_dir: Output directory
            target_bed: Target regions BED file
            bait_intervals: Bait interval list (Picard format)
            target_intervals: Target interval list (Picard format)

        Returns:
            Dictionary of pypipegraph jobs
        """
        output_dir = Path(output_dir)
        sample_dir = output_dir / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)

        jobs = {}

        # Alignment metrics
        align_metrics = sample_dir / f"{sample_name}.alignment_metrics.txt"
        jobs["alignment"] = self.collect_alignment_metrics(
            bam_file, reference_fasta, align_metrics
        )

        # Insert size metrics
        insert_metrics = sample_dir / f"{sample_name}.insert_size_metrics.txt"
        insert_histogram = sample_dir / f"{sample_name}.insert_size_histogram.pdf"
        jobs["insert_size"] = self.collect_insert_size_metrics(
            bam_file, insert_metrics, insert_histogram
        )

        # Coverage with mosdepth
        if target_bed:
            mosdepth_prefix = sample_dir / f"{sample_name}"
            jobs["coverage"] = self.run_mosdepth(bam_file, mosdepth_prefix, target_bed)

        # HS metrics (if intervals provided)
        if bait_intervals and target_intervals:
            hs_metrics = sample_dir / f"{sample_name}.hs_metrics.txt"
            jobs["hs_metrics"] = self.collect_hs_metrics(
                bam_file, reference_fasta, bait_intervals, target_intervals, hs_metrics
            )

        return jobs
