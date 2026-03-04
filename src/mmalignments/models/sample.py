# """Sample data models."""

# from dataclasses import dataclass, field
# from pathlib import Path
# from typing import Optional, List


# @dataclass
# class ExomeSeqSample:
#     """Represents an exome sequencing sample."""

#     sample_name: str
#     sample_type: str  # "tumor" or "normal"
#     fastq_r1: Optional[Path] = None
#     fastq_r2: Optional[Path] = None
#     bam_file: Optional[Path] = None
#     dedup_bam: Optional[Path] = None

#     # Metadata
#     mouse_id: Optional[str] = None
#     tissue: Optional[str] = None  # e.g., "Kidney", "SCC"
#     library_prep: Optional[str] = None
#     sequencing_platform: Optional[str] = None

#     # Paired sample
#     paired_sample: Optional[str] = None  # Name of matched tumor/normal


# @dataclass
# class TumorNormalPair:
#     """Represents a tumor-normal pair for somatic variant calling."""

#     tumor_sample: str
#     normal_sample: str
#     pair_id: Optional[str] = None

#     # Input files
#     tumor_bam: Optional[Path] = None
#     normal_bam: Optional[Path] = None

#     # Output files
#     unfiltered_vcf: Optional[Path] = None
#     filtered_vcf: Optional[Path] = None
#     annotated_vcf: Optional[Path] = None

#     # Intermediate files
#     f1r2_tar: Optional[Path] = None  # LearnReadOrientationModel
#     orientation_model: Optional[Path] = None
#     tumor_pileup: Optional[Path] = None
#     normal_pileup: Optional[Path] = None
#     contamination_table: Optional[Path] = None
#     segments_table: Optional[Path] = None

#     @property
#     def name(self) -> str:
#         """Get pair identifier."""
#         if self.pair_id:
#             return self.pair_id
#         return f"{self.tumor_sample}_vs_{self.normal_sample}"
