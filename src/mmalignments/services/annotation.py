"""Variant annotation services using VEP and VCF utilities."""

import pypipegraph2 as ppg
from pathlib import Path
from typing import Optional, List, Dict
import subprocess
import gzip
from ..models.variant import Variant, VariantType, VariantEffect


class AnnotationService:
    """Service for variant annotation."""

    def __init__(self, threads: int = 4):
        """
        Initialize annotation service.
        
        Args:
            threads: Number of threads for VEP
        """
        self.threads = threads

    def run_vep(
        self,
        input_vcf: Path,
        output_vcf: Path,
        reference_fasta: Path,
        cache_dir: Optional[Path] = None,
        species: str = "mus_musculus",
        assembly: str = "GRCm39",
        plugins: Optional[List[str]] = None,
    ) -> ppg.FileGeneratingJob:
        """
        Run Ensembl Variant Effect Predictor (VEP).
        
        Args:
            input_vcf: Input VCF file
            output_vcf: Output annotated VCF
            reference_fasta: Reference genome FASTA
            cache_dir: VEP cache directory
            species: Species name
            assembly: Genome assembly
            plugins: VEP plugins to use (e.g., ['CADD', 'REVEL'])
            
        Returns:
            pypipegraph job
        """
        output_vcf = Path(output_vcf)
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        
        def run_vep_cmd():
            cmd = [
                "vep",
                "--input_file", str(input_vcf),
                "--output_file", str(output_vcf),
                "--fasta", str(reference_fasta),
                "--species", species,
                "--assembly", assembly,
                "--vcf",
                "--force_overwrite",
                "--fork", str(self.threads),
                "--everything",  # Include all annotations
            ]
            
            if cache_dir:
                cmd.extend(["--dir_cache", str(cache_dir), "--cache"])
            else:
                cmd.append("--database")  # Use database instead of cache
            
            if plugins:
                for plugin in plugins:
                    cmd.extend(["--plugin", plugin])
            
            subprocess.run(cmd, check=True)

        job = ppg.FileGeneratingJob(output_vcf, run_vep_cmd)
        return job

    def filter_vcf_by_bed(
        self,
        input_vcf: Path,
        output_vcf: Path,
        bed_file: Path,
    ) -> ppg.FileGeneratingJob:
        """
        Filter VCF to only variants within BED regions using bcftools.
        
        Args:
            input_vcf: Input VCF file
            output_vcf: Output filtered VCF
            bed_file: BED file with target regions
            
        Returns:
            pypipegraph job
        """
        output_vcf = Path(output_vcf)
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        
        def filter_vcf():
            cmd = [
                "bcftools", "view",
                "-R", str(bed_file),
                "-O", "z",  # compressed output
                "-o", str(output_vcf),
                str(input_vcf),
            ]
            subprocess.run(cmd, check=True)
            
            # Index output VCF
            subprocess.run(["bcftools", "index", "-t", str(output_vcf)], check=True)

        job = ppg.FileGeneratingJob(output_vcf, filter_vcf)
        return job

    def filter_vcf_by_quality(
        self,
        input_vcf: Path,
        output_vcf: Path,
        min_dp: int = 10,
        min_alt_reads: int = 3,
        min_vaf: float = 0.0,
        pass_only: bool = True,
    ) -> ppg.FileGeneratingJob:
        """
        Filter VCF by quality thresholds using bcftools.
        
        Args:
            input_vcf: Input VCF file
            output_vcf: Output filtered VCF
            min_dp: Minimum depth
            min_alt_reads: Minimum alternate allele reads
            min_vaf: Minimum variant allele frequency
            pass_only: Keep only PASS variants
            
        Returns:
            pypipegraph job
        """
        output_vcf = Path(output_vcf)
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        
        def filter_vcf():
            # Build filter expression
            filters = []
            if pass_only:
                filters.append('FILTER="PASS"')
            if min_dp > 0:
                filters.append(f'FORMAT/DP>={min_dp}')
            if min_alt_reads > 0:
                filters.append(f'FORMAT/AD[0:1]>={min_alt_reads}')
            if min_vaf > 0:
                filters.append(f'(FORMAT/AD[0:1]/FORMAT/DP)>={min_vaf}')
            
            filter_expr = " && ".join(filters)
            
            cmd = [
                "bcftools", "view",
                "-i", filter_expr,
                "-O", "z",
                "-o", str(output_vcf),
                str(input_vcf),
            ]
            subprocess.run(cmd, check=True)
            
            # Index output VCF
            subprocess.run(["bcftools", "index", "-t", str(output_vcf)], check=True)

        job = ppg.FileGeneratingJob(output_vcf, filter_vcf)
        return job

    def parse_vcf_variants(
        self,
        vcf_file: Path,
        sample_name: Optional[str] = None,
    ) -> List[Variant]:
        """
        Parse VCF file and extract variants.
        
        Args:
            vcf_file: VCF file path
            sample_name: Sample name to extract (if None, use first sample)
            
        Returns:
            List of Variant objects
        """
        variants = []
        
        open_func = gzip.open if str(vcf_file).endswith('.gz') else open
        
        with open_func(vcf_file, 'rt') as f:
            samples = []
            for line in f:
                if line.startswith('##'):
                    continue
                
                if line.startswith('#CHROM'):
                    # Header line
                    parts = line.strip().split('\t')
                    samples = parts[9:]  # Sample names
                    if sample_name and sample_name in samples:
                        sample_idx = samples.index(sample_name)
                    else:
                        sample_idx = 0  # First sample
                    continue
                
                # Variant line
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                filter_status = parts[6]
                info = parts[7]
                
                # Determine variant type
                if len(ref) == 1 and len(alt) == 1:
                    variant_type = VariantType.SNV
                elif len(ref) < len(alt):
                    variant_type = VariantType.INSERTION
                elif len(ref) > len(alt):
                    variant_type = VariantType.DELETION
                else:
                    variant_type = VariantType.COMPLEX
                
                # Parse FORMAT and sample data
                format_fields = parts[8].split(':') if len(parts) > 8 else []
                sample_data = parts[9 + sample_idx].split(':') if len(parts) > 9 + sample_idx else []
                
                # Extract genotype information
                tumor_dp = 0
                tumor_ad_ref = 0
                tumor_ad_alt = 0
                
                if 'DP' in format_fields and len(sample_data) > format_fields.index('DP'):
                    try:
                        tumor_dp = int(sample_data[format_fields.index('DP')])
                    except (ValueError, IndexError):
                        pass
                
                if 'AD' in format_fields and len(sample_data) > format_fields.index('AD'):
                    try:
                        ad_str = sample_data[format_fields.index('AD')]
                        ad_parts = ad_str.split(',')
                        if len(ad_parts) >= 2:
                            tumor_ad_ref = int(ad_parts[0])
                            tumor_ad_alt = int(ad_parts[1])
                    except (ValueError, IndexError):
                        pass
                
                # Create variant object
                variant = Variant(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    variant_type=variant_type,
                    filter_status=filter_status,
                    tumor_dp=tumor_dp,
                    tumor_ad_ref=tumor_ad_ref,
                    tumor_ad_alt=tumor_ad_alt,
                )
                
                variants.append(variant)
        
        return variants

    def count_variants_by_effect(
        self,
        annotated_vcf: Path,
    ) -> Dict[str, int]:
        """
        Count variants by predicted effect from VEP-annotated VCF.
        
        Args:
            annotated_vcf: VEP-annotated VCF file
            
        Returns:
            Dictionary with counts per effect type
        """
        counts = {
            "total": 0,
            "coding": 0,
            "nonsynonymous": 0,
            "synonymous": 0,
            "missense": 0,
            "nonsense": 0,
            "frameshift": 0,
            "splice": 0,
            "intronic": 0,
            "intergenic": 0,
        }
        
        open_func = gzip.open if str(annotated_vcf).endswith('.gz') else open
        
        with open_func(annotated_vcf, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                counts["total"] += 1
                
                # Parse CSQ field from VEP
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                info = parts[7]
                csq_match = None
                for field in info.split(';'):
                    if field.startswith('CSQ='):
                        csq_match = field[4:]
                        break
                
                if not csq_match:
                    continue
                
                # Parse first consequence
                csq_parts = csq_match.split(',')[0].split('|')
                if len(csq_parts) > 1:
                    consequence = csq_parts[1].lower()
                    
                    if 'missense' in consequence:
                        counts["missense"] += 1
                        counts["nonsynonymous"] += 1
                        counts["coding"] += 1
                    elif 'nonsense' in consequence or 'stop_gained' in consequence:
                        counts["nonsense"] += 1
                        counts["nonsynonymous"] += 1
                        counts["coding"] += 1
                    elif 'frameshift' in consequence:
                        counts["frameshift"] += 1
                        counts["nonsynonymous"] += 1
                        counts["coding"] += 1
                    elif 'synonymous' in consequence:
                        counts["synonymous"] += 1
                        counts["coding"] += 1
                    elif 'splice' in consequence:
                        counts["splice"] += 1
                    elif 'intron' in consequence:
                        counts["intronic"] += 1
                    elif 'intergenic' in consequence:
                        counts["intergenic"] += 1
        
        return counts
