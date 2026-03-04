"""Callable regions calculation for tumor mutational burden."""

import pypipegraph2 as ppg
from pathlib import Path
from typing import Optional
import subprocess
import pysam


class CallableRegions:
    """Calculate callable bases for mutational burden analysis."""

    def __init__(self, min_depth: int = 10, min_base_quality: int = 20, min_mapping_quality: int = 20):
        """
        Initialize callable regions calculator.
        
        Args:
            min_depth: Minimum depth for callable bases
            min_base_quality: Minimum base quality
            min_mapping_quality: Minimum mapping quality
        """
        self.min_depth = min_depth
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality

    def calculate_callable_loci_gatk(
        self,
        bam_file: Path,
        reference_fasta: Path,
        output_bed: Path,
        intervals: Optional[Path] = None,
    ) -> ppg.FileGeneratingJob:
        """
        Calculate callable loci using GATK CallableLoci (deprecated but functional).
        
        Args:
            bam_file: Input BAM file
            reference_fasta: Reference genome
            output_bed: Output BED file with callable regions
            intervals: Target intervals
            
        Returns:
            pypipegraph job
        """
        output_bed = Path(output_bed)
        output_bed.parent.mkdir(parents=True, exist_ok=True)
        
        summary_file = output_bed.parent / f"{output_bed.stem}_summary.txt"
        
        def calculate():
            cmd = [
                "gatk", "CallableLoci",
                "-R", str(reference_fasta),
                "-I", str(bam_file),
                "--summary", str(summary_file),
                "-O", str(output_bed),
                f"--min-depth", str(self.min_depth),
                f"--min-base-quality", str(self.min_base_quality),
                f"--min-mapping-quality", str(self.min_mapping_quality),
            ]
            
            if intervals:
                cmd.extend(["-L", str(intervals)])
            
            subprocess.run(cmd, check=True)

        job = ppg.FileGeneratingJob(output_bed, calculate)
        return job

    def calculate_callable_from_mosdepth(
        self,
        mosdepth_summary: Path,
        mosdepth_regions: Optional[Path],
        output_file: Path,
    ) -> ppg.FileGeneratingJob:
        """
        Calculate callable bases from mosdepth output.
        
        Args:
            mosdepth_summary: mosdepth summary file
            mosdepth_regions: mosdepth regions file (optional)
            output_file: Output summary file
            
        Returns:
            pypipegraph job
        """
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        def calculate():
            total_bases = 0
            callable_bases = 0
            
            # Parse mosdepth summary
            with open(mosdepth_summary, 'r') as f:
                for line in f:
                    if line.startswith('total') or line.startswith('chr'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            region = parts[0]
                            length = int(parts[1])
                            bases_covered = int(parts[2])
                            mean_depth = float(parts[3])
                            
                            total_bases += length
                            
                            # Estimate callable bases (depth >= min_depth)
                            # This is approximate; for exact counts, parse regions file
                            if mean_depth >= self.min_depth:
                                callable_bases += bases_covered
            
            # If regions file provided, get more accurate counts
            if mosdepth_regions and Path(mosdepth_regions).exists():
                callable_bases = self._count_callable_from_regions(mosdepth_regions)
            
            # Write summary
            with open(output_file, 'w') as out:
                out.write(f"total_bases\t{total_bases}\n")
                out.write(f"callable_bases\t{callable_bases}\n")
                out.write(f"callable_mb\t{callable_bases / 1e6:.2f}\n")
                out.write(f"pct_callable\t{100 * callable_bases / total_bases if total_bases > 0 else 0:.2f}\n")

        job = ppg.FileGeneratingJob(output_file, calculate)
        return job

    def _count_callable_from_regions(self, regions_file: Path) -> int:
        """Count callable bases from mosdepth regions file."""
        callable_bases = 0
        
        # mosdepth regions file is BED format: chr start end depth
        import gzip
        open_func = gzip.open if str(regions_file).endswith('.gz') else open
        
        with open_func(regions_file, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    start = int(parts[1])
                    end = int(parts[2])
                    depth = float(parts[3])
                    
                    if depth >= self.min_depth:
                        callable_bases += (end - start)
        
        return callable_bases

    def count_callable_bases_in_bed(
        self,
        bam_file: Path,
        bed_file: Path,
        output_file: Path,
    ) -> ppg.FileGeneratingJob:
        """
        Count callable bases within BED regions using pysam.
        
        Args:
            bam_file: Input BAM file
            bed_file: Target BED file
            output_file: Output summary file
            
        Returns:
            pypipegraph job
        """
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        def count():
            total_bases = 0
            callable_bases = 0
            
            # Read BED file
            regions = []
            with open(bed_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.startswith('track'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        regions.append((chrom, start, end))
                        total_bases += (end - start)
            
            # Open BAM and count callable bases
            bamfile = pysam.AlignmentFile(str(bam_file), 'rb')
            
            for chrom, start, end in regions:
                for pileup_column in bamfile.pileup(
                    chrom, start, end,
                    truncate=True,
                    min_base_quality=self.min_base_quality,
                    min_mapping_quality=self.min_mapping_quality,
                ):
                    if pileup_column.n >= self.min_depth:
                        callable_bases += 1
            
            bamfile.close()
            
            # Write summary
            with open(output_file, 'w') as out:
                out.write(f"total_bases\t{total_bases}\n")
                out.write(f"callable_bases\t{callable_bases}\n")
                out.write(f"callable_mb\t{callable_bases / 1e6:.2f}\n")
                out.write(f"pct_callable\t{100 * callable_bases / total_bases if total_bases > 0 else 0:.2f}\n")

        job = ppg.FileGeneratingJob(output_file, count)
        return job

    def read_callable_summary(self, summary_file: Path) -> dict:
        """
        Read callable bases summary file.
        
        Args:
            summary_file: Summary file from callable calculation
            
        Returns:
            Dictionary with callable bases statistics
        """
        result = {}
        with open(summary_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    key = parts[0]
                    try:
                        value = float(parts[1])
                        result[key] = value
                    except ValueError:
                        result[key] = parts[1]
        return result
