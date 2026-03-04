"""BAM preprocessing: sorting, marking duplicates, indexing, and BQSR."""

import pypipegraph2 as ppg
from pathlib import Path
from typing import Dict, Optional, List
import subprocess


class BamPreprocessor:
    """Handles BAM file preprocessing steps."""

    def __init__(self, java_mem: str = "8g", threads: int = 4):
        """
        Initialize BAM preprocessor.
        
        Args:
            java_mem: Memory allocation for Java tools (e.g., "8g")
            threads: Number of threads to use
        """
        self.java_mem = java_mem
        self.threads = threads

    def sort_bam(
        self,
        input_bam: Path,
        output_bam: Path,
        temp_dir: Optional[Path] = None,
    ) -> ppg.FileGeneratingJob:
        """
        Sort BAM file by coordinate using samtools.
        
        Args:
            input_bam: Path to input BAM
            output_bam: Path to output sorted BAM
            temp_dir: Temporary directory for sorting
            
        Returns:
            pypipegraph job
        """
        output_bam = Path(output_bam)
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        
        def sort():
            cmd = [
                "samtools", "sort",
                "-@", str(self.threads),
                "-o", str(output_bam),
            ]
            if temp_dir:
                temp_dir.mkdir(parents=True, exist_ok=True)
                cmd.extend(["-T", str(temp_dir / "sort_tmp")])
            cmd.append(str(input_bam))
            
            subprocess.run(cmd, check=True)

        job = ppg.FileGeneratingJob(output_bam, sort)
        return job

    def index_bam(self, bam_file: Path) -> ppg.FileGeneratingJob:
        """
        Index BAM file using samtools.
        
        Args:
            bam_file: Path to BAM file to index
            
        Returns:
            pypipegraph job creating .bai file
        """
        bam_file = Path(bam_file)
        index_file = Path(str(bam_file) + ".bai")
        
        def index():
            cmd = [
                "samtools", "index",
                "-@", str(self.threads),
                str(bam_file),
            ]
            subprocess.run(cmd, check=True)

        job = ppg.FileGeneratingJob(index_file, index)
        return job

    def mark_duplicates(
        self,
        input_bam: Path,
        output_bam: Path,
        metrics_file: Path,
        remove_duplicates: bool = False,
        optical_duplicate_pixel_distance: int = 100,
    ) -> ppg.MultiFileGeneratingJob:
        """
        Mark duplicate reads using Picard MarkDuplicates.
        
        Args:
            input_bam: Path to input BAM
            output_bam: Path to output deduplicated BAM
            metrics_file: Path to metrics output file
            remove_duplicates: If True, remove duplicates instead of marking
            optical_duplicate_pixel_distance: Distance for optical duplicates
            
        Returns:
            pypipegraph job creating output BAM and metrics
        """
        output_bam = Path(output_bam)
        metrics_file = Path(metrics_file)
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        metrics_file.parent.mkdir(parents=True, exist_ok=True)
        
        def mark_dups():
            cmd = [
                "picard",
                f"-Xmx{self.java_mem}",
                "MarkDuplicates",
                f"INPUT={input_bam}",
                f"OUTPUT={output_bam}",
                f"METRICS_FILE={metrics_file}",
                f"REMOVE_DUPLICATES={'true' if remove_duplicates else 'false'}",
                "ASSUME_SORT_ORDER=coordinate",
                f"OPTICAL_DUPLICATE_PIXEL_DISTANCE={optical_duplicate_pixel_distance}",
                "CREATE_INDEX=true",
                "VALIDATION_STRINGENCY=LENIENT",
            ]
            subprocess.run(cmd, check=True)

        job = ppg.MultiFileGeneratingJob(
            [output_bam, metrics_file],
            mark_dups,
        )
        return job

    def base_quality_recalibration(
        self,
        input_bam: Path,
        output_bam: Path,
        reference_fasta: Path,
        known_sites_vcf: Optional[List[Path]] = None,
        intervals: Optional[Path] = None,
    ) -> ppg.MultiFileGeneratingJob:
        """
        Perform Base Quality Score Recalibration (BQSR) using GATK.
        
        Args:
            input_bam: Path to input BAM
            output_bam: Path to output recalibrated BAM
            reference_fasta: Path to reference genome FASTA
            known_sites_vcf: List of known variant VCF files (e.g., dbSNP)
            intervals: BED file for targeted regions
            
        Returns:
            pypipegraph job
        """
        output_bam = Path(output_bam)
        recal_table = output_bam.parent / f"{output_bam.stem}_recal_data.table"
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        
        def bqsr():
            # Step 1: BaseRecalibrator
            cmd_recal = [
                "gatk", "BaseRecalibrator",
                "-R", str(reference_fasta),
                "-I", str(input_bam),
                "-O", str(recal_table),
            ]
            
            if known_sites_vcf:
                for vcf in known_sites_vcf:
                    cmd_recal.extend(["--known-sites", str(vcf)])
            
            if intervals:
                cmd_recal.extend(["-L", str(intervals)])
            
            subprocess.run(cmd_recal, check=True)
            
            # Step 2: ApplyBQSR
            cmd_apply = [
                "gatk", "ApplyBQSR",
                "-R", str(reference_fasta),
                "-I", str(input_bam),
                "--bqsr-recal-file", str(recal_table),
                "-O", str(output_bam),
            ]
            
            subprocess.run(cmd_apply, check=True)

        job = ppg.MultiFileGeneratingJob(
            [output_bam, recal_table],
            bqsr,
        )
        return job

    def preprocess_bam(
        self,
        input_bam: Path,
        output_dir: Path,
        sample_name: str,
        reference_fasta: Optional[Path] = None,
        known_sites_vcf: Optional[List[Path]] = None,
        intervals: Optional[Path] = None,
        skip_bqsr: bool = True,  # BQSR requires known sites, skip by default for mouse
    ) -> Dict[str, ppg.Job]:
        """
        Complete BAM preprocessing pipeline.
        
        Steps:
        1. Sort by coordinate
        2. Mark duplicates
        3. Index
        4. (Optional) BQSR
        5. Final index
        
        Args:
            input_bam: Path to input BAM
            output_dir: Output directory
            sample_name: Sample name
            reference_fasta: Reference genome (required for BQSR)
            known_sites_vcf: Known variant sites (required for BQSR)
            intervals: Target intervals
            skip_bqsr: Skip BQSR step (default True for mouse)
            
        Returns:
            Dictionary of jobs with keys: 'sort', 'markdup', 'index', 'bqsr', 'final_bam'
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        jobs = {}
        
        # Sort
        sorted_bam = output_dir / f"{sample_name}.sorted.bam"
        jobs['sort'] = self.sort_bam(input_bam, sorted_bam)
        
        # Mark duplicates
        dedup_bam = output_dir / f"{sample_name}.dedup.bam"
        metrics_file = output_dir / f"{sample_name}.dedup.metrics.txt"
        jobs['markdup'] = self.mark_duplicates(sorted_bam, dedup_bam, metrics_file)
        jobs['markdup'].depends_on(jobs['sort'])
        
        # Index deduplicated BAM
        jobs['index_dedup'] = self.index_bam(dedup_bam)
        jobs['index_dedup'].depends_on(jobs['markdup'])
        
        # BQSR (optional)
        if not skip_bqsr and reference_fasta:
            recal_bam = output_dir / f"{sample_name}.dedup.recal.bam"
            jobs['bqsr'] = self.base_quality_recalibration(
                dedup_bam, recal_bam, reference_fasta, known_sites_vcf, intervals
            )
            jobs['bqsr'].depends_on(jobs['index_dedup'])
            
            # Index final BAM
            jobs['index_final'] = self.index_bam(recal_bam)
            jobs['index_final'].depends_on(jobs['bqsr'])
            
            jobs['final_bam'] = recal_bam
        else:
            jobs['final_bam'] = dedup_bam
        
        return jobs
