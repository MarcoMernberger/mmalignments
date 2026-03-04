"""Utility functions for BED file operations."""

import pypipegraph2 as ppg
from pathlib import Path
import subprocess
from typing import Optional


def create_padded_bed(
    input_bed: Path,
    output_bed: Path,
    padding: int = 100,
    genome_file: Optional[Path] = None,
) -> ppg.FileGeneratingJob:
    """
    Create padded BED file using bedtools slop.
    
    Args:
        input_bed: Input BED file
        output_bed: Output padded BED file
        padding: Padding in bases on each side
        genome_file: Genome file with chromosome sizes (required)
        
    Returns:
        pypipegraph job
    """
    output_bed = Path(output_bed)
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    
    def create_padded():
        if genome_file is None:
            raise ValueError("genome_file is required for creating padded BED")
        
        cmd = [
            "bedtools", "slop",
            "-i", str(input_bed),
            "-g", str(genome_file),
            "-b", str(padding),
        ]
        
        with open(output_bed, 'w') as out:
            subprocess.run(cmd, stdout=out, check=True)

    job = ppg.FileGeneratingJob(output_bed, create_padded)
    return job


def bed_to_interval_list(
    bed_file: Path,
    reference_dict: Path,
    output_interval_list: Path,
) -> ppg.FileGeneratingJob:
    """
    Convert BED file to Picard interval list format.
    
    Args:
        bed_file: Input BED file
        reference_dict: Reference sequence dictionary (.dict file)
        output_interval_list: Output interval list
        
    Returns:
        pypipegraph job
    """
    output_interval_list = Path(output_interval_list)
    output_interval_list.parent.mkdir(parents=True, exist_ok=True)
    
    def convert():
        cmd = [
            "picard", "BedToIntervalList",
            f"INPUT={bed_file}",
            f"SEQUENCE_DICTIONARY={reference_dict}",
            f"OUTPUT={output_interval_list}",
        ]
        subprocess.run(cmd, check=True)

    job = ppg.FileGeneratingJob(output_interval_list, convert)
    return job


def create_genome_file_from_fasta(
    fasta_file: Path,
    output_genome_file: Path,
) -> ppg.FileGeneratingJob:
    """
    Create genome file (chromosome sizes) from FASTA using samtools.
    
    Args:
        fasta_file: Input FASTA file
        output_genome_file: Output genome file
        
    Returns:
        pypipegraph job
    """
    output_genome_file = Path(output_genome_file)
    output_genome_file.parent.mkdir(parents=True, exist_ok=True)
    
    def create_genome():
        # Get fai index
        fai_file = Path(str(fasta_file) + ".fai")
        if not fai_file.exists():
            subprocess.run(["samtools", "faidx", str(fasta_file)], check=True)
        
        # Extract chromosome sizes from .fai
        with open(fai_file, 'r') as f_in, open(output_genome_file, 'w') as f_out:
            for line in f_in:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom = parts[0]
                    size = parts[1]
                    f_out.write(f"{chrom}\t{size}\n")

    job = ppg.FileGeneratingJob(output_genome_file, create_genome)
    return job


def merge_bed_files(
    bed_files: list[Path],
    output_bed: Path,
    sorted: bool = True,
) -> ppg.FileGeneratingJob:
    """
    Merge multiple BED files using bedtools merge.
    
    Args:
        bed_files: List of input BED files
        output_bed: Output merged BED file
        sorted: Whether to sort before merging
        
    Returns:
        pypipegraph job
    """
    output_bed = Path(output_bed)
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    
    def merge():
        if sorted:
            # Concatenate and sort
            temp_concat = output_bed.parent / f"{output_bed.stem}_temp.bed"
            with open(temp_concat, 'w') as out:
                for bed in bed_files:
                    with open(bed, 'r') as f:
                        out.write(f.read())
            
            # Sort
            temp_sorted = output_bed.parent / f"{output_bed.stem}_temp_sorted.bed"
            cmd_sort = ["bedtools", "sort", "-i", str(temp_concat)]
            with open(temp_sorted, 'w') as out:
                subprocess.run(cmd_sort, stdout=out, check=True)
            
            # Merge
            cmd_merge = ["bedtools", "merge", "-i", str(temp_sorted)]
            with open(output_bed, 'w') as out:
                subprocess.run(cmd_merge, stdout=out, check=True)
            
            # Cleanup
            temp_concat.unlink()
            temp_sorted.unlink()
        else:
            # Direct merge
            cmd = ["bedtools", "merge", "-i"]
            cmd.extend([str(b) for b in bed_files])
            with open(output_bed, 'w') as out:
                subprocess.run(cmd, stdout=out, check=True)

    job = ppg.FileGeneratingJob(output_bed, merge)
    return job
