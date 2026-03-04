# """Mutation and variant data models."""

# from dataclasses import dataclass, field
# from enum import Enum
# from pathlib import Path
# from typing import Any, Dict, List, Optional


# class VariantType(Enum):
#     """Type of genetic variant."""

#     SNV = "SNV"
#     INSERTION = "INSERTION"
#     DELETION = "DELETION"
#     MNV = "MNV"  # Multi-nucleotide variant
#     COMPLEX = "COMPLEX"


# class VariantEffect(Enum):
#     """Predicted effect of variant."""

#     SYNONYMOUS = "synonymous"
#     MISSENSE = "missense"
#     NONSENSE = "nonsense"
#     FRAMESHIFT = "frameshift"
#     INFRAME_INSERTION = "inframe_insertion"
#     INFRAME_DELETION = "inframe_deletion"
#     SPLICE_SITE = "splice_site"
#     SPLICE_REGION = "splice_region"
#     UTR_5 = "5_prime_UTR"
#     UTR_3 = "3_prime_UTR"
#     INTRON = "intron"
#     INTERGENIC = "intergenic"
#     UNKNOWN = "unknown"


# @dataclass
# class Variant:
#     """Represents a single genetic variant."""

#     chrom: str
#     pos: int
#     ref: str
#     alt: str
#     variant_type: VariantType
#     filter_status: str  # PASS, germline, weak_evidence, etc.

#     # Genotype information
#     tumor_dp: int = 0  # Tumor depth
#     tumor_ad_ref: int = 0  # Tumor allelic depth REF
#     tumor_ad_alt: int = 0  # Tumor allelic depth ALT
#     normal_dp: int = 0  # Normal depth
#     normal_ad_ref: int = 0
#     normal_ad_alt: int = 0

#     # Annotations
#     gene: Optional[str] = None
#     transcript: Optional[str] = None
#     effect: Optional[VariantEffect] = None
#     protein_change: Optional[str] = None
#     cdna_change: Optional[str] = None

#     # Quality metrics
#     quality_score: Optional[float] = None
#     strand_bias: Optional[float] = None
#     mapping_quality: Optional[float] = None

#     # Population frequency (if annotated)
#     gnomad_af: Optional[float] = None
#     gnomad_af_popmax: Optional[float] = None

#     # Prediction scores
#     cadd_score: Optional[float] = None
#     revel_score: Optional[float] = None
#     spliceai_score: Optional[float] = None

#     # Additional info
#     info: Dict[str, Any] = field(default_factory=dict)

#     @property
#     def tumor_vaf(self) -> float:
#         """Calculate tumor variant allele frequency."""
#         if self.tumor_dp == 0:
#             return 0.0
#         return self.tumor_ad_alt / self.tumor_dp

#     @property
#     def normal_vaf(self) -> float:
#         """Calculate normal variant allele frequency."""
#         if self.normal_dp == 0:
#             return 0.0
#         return self.normal_ad_alt / self.normal_dp

#     @property
#     def is_pass(self) -> bool:
#         """Check if variant passed all filters."""
#         return self.filter_status == "PASS"

#     @property
#     def is_coding(self) -> bool:
#         """Check if variant is in coding region."""
#         if not self.effect:
#             return False
#         coding_effects = {
#             VariantEffect.SYNONYMOUS,
#             VariantEffect.MISSENSE,
#             VariantEffect.NONSENSE,
#             VariantEffect.FRAMESHIFT,
#             VariantEffect.INFRAME_INSERTION,
#             VariantEffect.INFRAME_DELETION,
#         }
#         return self.effect in coding_effects

#     @property
#     def is_nonsynonymous(self) -> bool:
#         """Check if variant changes protein sequence."""
#         if not self.effect:
#             return False
#         nonsynonymous_effects = {
#             VariantEffect.MISSENSE,
#             VariantEffect.NONSENSE,
#             VariantEffect.FRAMESHIFT,
#             VariantEffect.INFRAME_INSERTION,
#             VariantEffect.INFRAME_DELETION,
#         }
#         return self.effect in nonsynonymous_effects

#     def to_dict(self) -> Dict[str, Any]:
#         """Convert variant to dictionary."""
#         return {
#             "chrom": self.chrom,
#             "pos": self.pos,
#             "ref": self.ref,
#             "alt": self.alt,
#             "variant_type": self.variant_type.value,
#             "filter": self.filter_status,
#             "tumor_dp": self.tumor_dp,
#             "tumor_vaf": self.tumor_vaf,
#             "normal_dp": self.normal_dp,
#             "normal_vaf": self.normal_vaf,
#             "gene": self.gene,
#             "effect": self.effect.value if self.effect else None,
#             "protein_change": self.protein_change,
#             "quality": self.quality_score,
#         }


# @dataclass
# class MutationalLoad:
#     """Mutational load (TMB) calculation for a tumor sample."""

#     sample_name: str
#     normal_name: str

#     # Variant counts
#     n_snv: int = 0
#     n_insertion: int = 0
#     n_deletion: int = 0
#     n_snv_coding: int = 0
#     n_snv_nonsynonymous: int = 0
#     n_indel_frameshift: int = 0
#     n_splice: int = 0

#     # Callable bases
#     callable_bases: int = 0
#     target_bases: int = 0

#     # Files
#     filtered_vcf: Optional[Path] = None
#     annotated_vcf: Optional[Path] = None

#     # QC metrics
#     mean_tumor_depth: Optional[float] = None
#     mean_normal_depth: Optional[float] = None
#     pct_callable: Optional[float] = None

#     @property
#     def n_total_variants(self) -> int:
#         """Total number of somatic variants."""
#         return self.n_snv + self.n_insertion + self.n_deletion

#     @property
#     def n_indels(self) -> int:
#         """Total number of indels."""
#         return self.n_insertion + self.n_deletion

#     @property
#     def callable_mb(self) -> float:
#         """Callable megabases."""
#         return self.callable_bases / 1e6

#     @property
#     def mutational_load(self) -> float:
#         """Mutational load (mutations per Mb)."""
#         if self.callable_mb == 0:
#             return 0.0
#         return self.n_total_variants / self.callable_mb

#     @property
#     def snv_load(self) -> float:
#         """SNV load (SNVs per Mb)."""
#         if self.callable_mb == 0:
#             return 0.0
#         return self.n_snv / self.callable_mb

#     @property
#     def indel_load(self) -> float:
#         """Indel load (indels per Mb)."""
#         if self.callable_mb == 0:
#             return 0.0
#         return self.n_indels / self.callable_mb

#     def to_dict(self) -> Dict[str, Any]:
#         """Convert to dictionary for reporting."""
#         return {
#             "sample": self.sample_name,
#             "normal": self.normal_name,
#             "n_snv": self.n_snv,
#             "n_insertion": self.n_insertion,
#             "n_deletion": self.n_deletion,
#             "n_total": self.n_total_variants,
#             "n_snv_coding": self.n_snv_coding,
#             "n_snv_nonsynonymous": self.n_snv_nonsynonymous,
#             "n_indel_frameshift": self.n_indel_frameshift,
#             "n_splice": self.n_splice,
#             "callable_bases": self.callable_bases,
#             "callable_mb": self.callable_mb,
#             "mutational_load": self.mutational_load,
#             "snv_load": self.snv_load,
#             "indel_load": self.indel_load,
#             "mean_tumor_depth": self.mean_tumor_depth,
#             "mean_normal_depth": self.mean_normal_depth,
#             "pct_callable": self.pct_callable,
#         }


# @dataclass
# class SampleMutationReport:
#     """Comprehensive mutation report for a tumor-normal pair."""

#     tumor_name: str
#     normal_name: str
#     mutational_load: MutationalLoad
#     variants: List[Variant] = field(default_factory=list)
#     qc_metrics: Optional[Dict[str, Any]] = None

#     def filter_variants(
#         self,
#         min_tumor_dp: int = 10,
#         min_alt_reads: int = 3,
#         min_vaf: float = 0.0,
#         pass_only: bool = True,
#     ) -> List[Variant]:
#         """Filter variants by quality thresholds."""
#         filtered = []
#         for var in self.variants:
#             if pass_only and not var.is_pass:
#                 continue
#             if var.tumor_dp < min_tumor_dp:
#                 continue
#             if var.tumor_ad_alt < min_alt_reads:
#                 continue
#             if var.tumor_vaf < min_vaf:
#                 continue
#             filtered.append(var)
#         return filtered

#     def get_coding_variants(self) -> List[Variant]:
#         """Get all coding variants."""
#         return [v for v in self.variants if v.is_coding]

#     def get_nonsynonymous_variants(self) -> List[Variant]:
#         """Get all nonsynonymous variants."""
#         return [v for v in self.variants if v.is_nonsynonymous]
