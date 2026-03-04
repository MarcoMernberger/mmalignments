"""Quality control metrics data models."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, Any


@dataclass
class FastQCMetrics:
    """FastQC quality metrics for a single FASTQ file."""

    sample_name: str
    read_type: str  # "R1" or "R2" or "single"
    total_sequences: int
    poor_quality: int
    sequence_length: str
    gc_content: int
    per_base_quality_passed: bool
    per_sequence_quality_passed: bool
    adapter_content_passed: bool
    duplication_level: Optional[float] = None
    overrepresented_sequences: int = 0
    fastqc_report_html: Optional[Path] = None
    fastqc_data_txt: Optional[Path] = None

    @property
    def quality_score(self) -> str:
        """Overall quality assessment."""
        if all([self.per_base_quality_passed, self.per_sequence_quality_passed, self.adapter_content_passed]):
            return "PASS"
        elif any([not self.per_base_quality_passed, not self.per_sequence_quality_passed]):
            return "FAIL"
        else:
            return "WARN"


@dataclass
class AlignmentMetrics:
    """Picard alignment summary metrics."""

    sample_name: str
    total_reads: int
    mapped_reads: int
    paired_reads: int
    properly_paired: int
    read_duplicates: int
    mean_coverage: float
    median_coverage: float
    pct_target_bases_10x: float
    pct_target_bases_20x: float
    pct_target_bases_30x: float
    fold_80_base_penalty: Optional[float] = None
    zero_cvg_targets_pct: Optional[float] = None
    pct_selected_bases: Optional[float] = None  # on-target %
    mean_target_coverage: Optional[float] = None
    median_target_coverage: Optional[float] = None

    @property
    def mapping_rate(self) -> float:
        """Calculate mapping rate."""
        if self.total_reads == 0:
            return 0.0
        return (self.mapped_reads / self.total_reads) * 100

    @property
    def duplication_rate(self) -> float:
        """Calculate duplication rate."""
        if self.mapped_reads == 0:
            return 0.0
        return (self.read_duplicates / self.mapped_reads) * 100

    @property
    def proper_pair_rate(self) -> float:
        """Calculate properly paired rate."""
        if self.paired_reads == 0:
            return 0.0
        return (self.properly_paired / self.paired_reads) * 100


@dataclass
class InsertSizeMetrics:
    """Picard insert size metrics for paired-end data."""

    sample_name: str
    median_insert_size: int
    mean_insert_size: float
    standard_deviation: float
    min_insert_size: int
    max_insert_size: int
    pair_orientation: str = "FR"  # FR, RF, TANDEM


@dataclass
class CoverageMetrics:
    """mosdepth coverage metrics on target regions."""

    sample_name: str
    target_bed: Path
    total_target_bases: int
    covered_bases: int
    mean_coverage: float
    median_coverage: float
    coverage_10x_bases: int
    coverage_20x_bases: int
    coverage_30x_bases: int
    mosdepth_summary: Optional[Path] = None
    mosdepth_regions: Optional[Path] = None

    @property
    def pct_covered(self) -> float:
        """Percentage of target bases with any coverage."""
        if self.total_target_bases == 0:
            return 0.0
        return (self.covered_bases / self.total_target_bases) * 100

    @property
    def pct_10x(self) -> float:
        """Percentage of target bases with >=10x coverage."""
        if self.total_target_bases == 0:
            return 0.0
        return (self.coverage_10x_bases / self.total_target_bases) * 100

    @property
    def pct_20x(self) -> float:
        """Percentage of target bases with >=20x coverage."""
        if self.total_target_bases == 0:
            return 0.0
        return (self.coverage_20x_bases / self.total_target_bases) * 100

    @property
    def pct_30x(self) -> float:
        """Percentage of target bases with >=30x coverage."""
        if self.total_target_bases == 0:
            return 0.0
        return (self.coverage_30x_bases / self.total_target_bases) * 100


@dataclass
class DuplicationMetrics:
    """Picard MarkDuplicates metrics."""

    sample_name: str
    unpaired_reads_examined: int
    read_pairs_examined: int
    unmapped_reads: int
    unpaired_read_duplicates: int
    read_pair_duplicates: int
    read_pair_optical_duplicates: int
    percent_duplication: float
    estimated_library_size: int
    metrics_file: Optional[Path] = None


@dataclass
class SampleQCReport:
    """Comprehensive QC report for a sample."""

    sample_name: str
    fastqc_r1: Optional[FastQCMetrics] = None
    fastqc_r2: Optional[FastQCMetrics] = None
    alignment: Optional[AlignmentMetrics] = None
    insert_size: Optional[InsertSizeMetrics] = None
    coverage: Optional[CoverageMetrics] = None
    duplication: Optional[DuplicationMetrics] = None
    additional_metrics: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for reporting."""
        result = {"sample_name": self.sample_name}

        if self.alignment:
            result.update({
                "mapping_rate": self.alignment.mapping_rate,
                "duplication_rate": self.alignment.duplication_rate,
                "mean_coverage": self.alignment.mean_coverage,
                "pct_target_10x": self.alignment.pct_target_bases_10x,
                "on_target_pct": self.alignment.pct_selected_bases or 0.0,
            })

        if self.insert_size:
            result.update({
                "median_insert_size": self.insert_size.median_insert_size,
                "mean_insert_size": self.insert_size.mean_insert_size,
            })

        if self.coverage:
            result.update({
                "mean_target_coverage": self.coverage.mean_coverage,
                "pct_10x": self.coverage.pct_10x,
                "pct_20x": self.coverage.pct_20x,
                "pct_30x": self.coverage.pct_30x,
            })

        if self.duplication:
            result.update({
                "percent_duplication": self.duplication.percent_duplication,
            })

        result.update(self.additional_metrics)

        return result
