"""GATK Mutect2 somatic variant calling pipeline."""

import pypipegraph2 as ppg
from pathlib import Path
from typing import Optional, List, Dict
import subprocess


class Mutect2Caller:
    """Handles GATK Mutect2 somatic variant calling."""

    def __init__(self, java_mem: str = "8g", threads: int = 4):
        """
        Initialize Mutect2 caller.
        
        Args:
            java_mem: Memory allocation for GATK (e.g., "8g")
            threads: Number of threads to use
        """
        self.java_mem = java_mem
        self.threads = threads

    def call_somatic_variants(
        self,
        tumor_bam: Path,
        normal_bam: Path,
        reference_fasta: Path,
        output_vcf: Path,
        tumor_sample: str,
        normal_sample: str,
        intervals: Optional[Path] = None,
        panel_of_normals: Optional[Path] = None,
        germline_resource: Optional[Path] = None,
        f1r2_tar: Optional[Path] = None,
    ) -> ppg.MultiFileGeneratingJob:
        """
        Run Mutect2 somatic variant calling.
        
        Args:
            tumor_bam: Path to tumor BAM file
            normal_bam: Path to normal BAM file
            reference_fasta: Path to reference genome
            output_vcf: Path to output VCF
            tumor_sample: Tumor sample name
            normal_sample: Normal sample name
            intervals: Target intervals BED file
            panel_of_normals: Panel of normals VCF
            germline_resource: Germline resource VCF (e.g., gnomAD)
            f1r2_tar: Output for F1R2 counts (for LearnReadOrientationModel)
            
        Returns:
            pypipegraph job
        """
        output_vcf = Path(output_vcf)
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        
        # Mutect2 also creates .vcf.gz.stats
        stats_file = Path(str(output_vcf) + ".stats")
        outputs = [output_vcf, stats_file]
        
        if f1r2_tar:
            f1r2_tar = Path(f1r2_tar)
            outputs.append(f1r2_tar)
        
        def call_variants():
            cmd = [
                "gatk", "Mutect2",
                "-R", str(reference_fasta),
                "-I", str(tumor_bam),
                "-I", str(normal_bam),
                "-tumor", tumor_sample,
                "-normal", normal_sample,
                "-O", str(output_vcf),
                "--native-pair-hmm-threads", str(self.threads),
            ]
            
            if intervals:
                cmd.extend(["-L", str(intervals)])
            
            if panel_of_normals:
                cmd.extend(["--panel-of-normals", str(panel_of_normals)])
            
            if germline_resource:
                cmd.extend(["--germline-resource", str(germline_resource)])
            
            if f1r2_tar:
                cmd.extend(["--f1r2-tar-gz", str(f1r2_tar)])
            
            subprocess.run(cmd, check=True)

        job = ppg.MultiFileGeneratingJob(outputs, call_variants)
        return job

    def learn_read_orientation_model(
        self,
        f1r2_tar: Path,
        output_model: Path,
    ) -> ppg.FileGeneratingJob:
        """
        Learn read orientation artifacts from F1R2 data.
        
        Args:
            f1r2_tar: F1R2 tar.gz from Mutect2
            output_model: Output orientation model
            
        Returns:
            pypipegraph job
        """
        output_model = Path(output_model)
        output_model.parent.mkdir(parents=True, exist_ok=True)
        
        def learn_model():
            cmd = [
                "gatk", "LearnReadOrientationModel",
                "-I", str(f1r2_tar),
                "-O", str(output_model),
            ]
            subprocess.run(cmd, check=True)

        job = ppg.FileGeneratingJob(output_model, learn_model)
        return job

    def get_pileup_summaries(
        self,
        bam_file: Path,
        reference_fasta: Path,
        output_table: Path,
        germline_resource: Path,
        intervals: Optional[Path] = None,
    ) -> ppg.FileGeneratingJob:
        """
        Get pileup summaries for contamination estimation.
        
        Args:
            bam_file: BAM file (tumor or normal)
            reference_fasta: Reference genome
            output_table: Output pileup table
            germline_resource: Common germline variants (e.g., gnomAD)
            intervals: Target intervals
            
        Returns:
            pypipegraph job
        """
        output_table = Path(output_table)
        output_table.parent.mkdir(parents=True, exist_ok=True)
        
        def get_pileups():
            cmd = [
                "gatk", "GetPileupSummaries",
                "-I", str(bam_file),
                "-R", str(reference_fasta),
                "-V", str(germline_resource),
                "-O", str(output_table),
            ]
            
            if intervals:
                cmd.extend(["-L", str(intervals)])
            
            subprocess.run(cmd, check=True)

        job = ppg.FileGeneratingJob(output_table, get_pileups)
        return job

    def calculate_contamination(
        self,
        tumor_pileup: Path,
        output_table: Path,
        normal_pileup: Optional[Path] = None,
        output_segments: Optional[Path] = None,
    ) -> ppg.MultiFileGeneratingJob:
        """
        Calculate contamination from pileup summaries.
        
        Args:
            tumor_pileup: Tumor pileup table
            output_table: Output contamination table
            normal_pileup: Normal pileup table (optional)
            output_segments: Output tumor segmentation table
            
        Returns:
            pypipegraph job
        """
        output_table = Path(output_table)
        output_table.parent.mkdir(parents=True, exist_ok=True)
        
        outputs = [output_table]
        if output_segments:
            output_segments = Path(output_segments)
            outputs.append(output_segments)
        
        def calc_contamination():
            cmd = [
                "gatk", "CalculateContamination",
                "-I", str(tumor_pileup),
                "-O", str(output_table),
            ]
            
            if normal_pileup:
                cmd.extend(["-matched", str(normal_pileup)])
            
            if output_segments:
                cmd.extend(["--tumor-segmentation", str(output_segments)])
            
            subprocess.run(cmd, check=True)

        job = ppg.MultiFileGeneratingJob(outputs, calc_contamination)
        return job

    def filter_mutect_calls(
        self,
        input_vcf: Path,
        reference_fasta: Path,
        output_vcf: Path,
        contamination_table: Optional[Path] = None,
        orientation_model: Optional[Path] = None,
        segments_table: Optional[Path] = None,
        stats_file: Optional[Path] = None,
    ) -> ppg.FileGeneratingJob:
        """
        Filter Mutect2 calls using contamination and orientation models.
        
        Args:
            input_vcf: Unfiltered VCF from Mutect2
            reference_fasta: Reference genome
            output_vcf: Filtered output VCF
            contamination_table: Contamination table
            orientation_model: Read orientation model
            segments_table: Tumor segmentation table
            stats_file: Mutect2 stats file
            
        Returns:
            pypipegraph job
        """
        output_vcf = Path(output_vcf)
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        
        def filter_calls():
            cmd = [
                "gatk", "FilterMutectCalls",
                "-R", str(reference_fasta),
                "-V", str(input_vcf),
                "-O", str(output_vcf),
            ]
            
            if contamination_table:
                cmd.extend(["--contamination-table", str(contamination_table)])
            
            if orientation_model:
                cmd.extend(["--ob-priors", str(orientation_model)])
            
            if segments_table:
                cmd.extend(["--tumor-segmentation", str(segments_table)])
            
            if stats_file:
                cmd.extend(["--stats", str(stats_file)])
            
            subprocess.run(cmd, check=True)

        job = ppg.FileGeneratingJob(output_vcf, filter_calls)
        return job

    def complete_mutect2_pipeline(
        self,
        tumor_bam: Path,
        normal_bam: Path,
        reference_fasta: Path,
        output_dir: Path,
        pair_name: str,
        tumor_sample: str,
        normal_sample: str,
        intervals: Optional[Path] = None,
        germline_resource: Optional[Path] = None,
        panel_of_normals: Optional[Path] = None,
    ) -> Dict[str, ppg.Job]:
        """
        Complete Mutect2 pipeline with all filtering steps.
        
        Args:
            tumor_bam: Tumor BAM file
            normal_bam: Normal BAM file
            reference_fasta: Reference genome
            output_dir: Output directory
            pair_name: Tumor-normal pair name
            tumor_sample: Tumor sample name
            normal_sample: Normal sample name
            intervals: Target intervals
            germline_resource: Germline resource VCF
            panel_of_normals: Panel of normals VCF
            
        Returns:
            Dictionary of pypipegraph jobs
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        jobs = {}
        
        # 1. Mutect2 calling
        unfiltered_vcf = output_dir / f"{pair_name}.unfiltered.vcf.gz"
        f1r2_tar = output_dir / f"{pair_name}.f1r2.tar.gz"
        
        jobs['mutect2'] = self.call_somatic_variants(
            tumor_bam, normal_bam, reference_fasta, unfiltered_vcf,
            tumor_sample, normal_sample, intervals, panel_of_normals,
            germline_resource, f1r2_tar
        )
        
        # 2. Learn read orientation model
        orientation_model = output_dir / f"{pair_name}.read_orientation_model.tar.gz"
        jobs['orientation'] = self.learn_read_orientation_model(f1r2_tar, orientation_model)
        jobs['orientation'].depends_on(jobs['mutect2'])
        
        # 3. Get pileup summaries (if germline resource available)
        if germline_resource:
            tumor_pileup = output_dir / f"{pair_name}.tumor_pileups.table"
            normal_pileup = output_dir / f"{pair_name}.normal_pileups.table"
            
            jobs['tumor_pileup'] = self.get_pileup_summaries(
                tumor_bam, reference_fasta, tumor_pileup, germline_resource, intervals
            )
            
            jobs['normal_pileup'] = self.get_pileup_summaries(
                normal_bam, reference_fasta, normal_pileup, germline_resource, intervals
            )
            
            # 4. Calculate contamination
            contamination_table = output_dir / f"{pair_name}.contamination.table"
            segments_table = output_dir / f"{pair_name}.segments.table"
            
            jobs['contamination'] = self.calculate_contamination(
                tumor_pileup, contamination_table, normal_pileup, segments_table
            )
            jobs['contamination'].depends_on(jobs['tumor_pileup'])
            jobs['contamination'].depends_on(jobs['normal_pileup'])
        else:
            contamination_table = None
            segments_table = None
        
        # 5. Filter Mutect calls
        filtered_vcf = output_dir / f"{pair_name}.filtered.vcf.gz"
        stats_file = Path(str(unfiltered_vcf) + ".stats")
        
        jobs['filter'] = self.filter_mutect_calls(
            unfiltered_vcf, reference_fasta, filtered_vcf,
            contamination_table, orientation_model, segments_table, stats_file
        )
        jobs['filter'].depends_on(jobs['mutect2'])
        jobs['filter'].depends_on(jobs['orientation'])
        if germline_resource:
            jobs['filter'].depends_on(jobs['contamination'])
        
        jobs['final_vcf'] = filtered_vcf
        
        return jobs
