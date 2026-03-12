"""Lightweight datamodels mirroring parts of MBF Sample and EnsemblGenome.

These are small, picklable dataclasses intended for use in mmalignments pipelines
where a compact representation of a sample or genome is sufficient.

They intentionally only include the most commonly used fields from
`mbf.align.raw.Sample` and the Ensembl genome wrapper in `mbf.genomes.ensembl`.
If you need more fields later, add them here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, List, Literal

# from .elements import Element


class Pairing(str, Enum):
    PAIRED = "paired"
    SINGLE = "single"


@dataclass
class Sample:

    name: str
    pairing: Literal["single", "paired"] = "single"
    reverse_reads: bool = False
    fastq_r1_path: Path | str | None = None
    fastq_r2_path: Path | str | None = None
    base_cache_dir: Path | str = Path("cache/lanes")
    base_result_dir: Path | str = Path("results/lanes")
    read_group: str | None = None
    cache_dir: Path = field(init=False)
    result_dir: Path = field(init=False)
    input_files: List[Path] = field(init=False)

    def __post_init__(self) -> None:
        self.base_cache_dir = Path(self.base_cache_dir)
        self.base_result_dir = Path(self.base_result_dir)

        self.cache_dir = self.base_cache_dir / self.name
        self.result_dir = self.base_result_dir / self.name

        # If user provided explicit FASTQs → use them
        if self.fastq_r1_path:
            r1 = Path(self.fastq_r1_path)
            r2 = Path(self.fastq_r2_path) if self.fastq_r2_path else None
        else:
            r1 = self.cache_dir / "input_R_1.fastq"
            r2 = (
                self.cache_dir / "input_R_2.fastq" if self.pairing == "paired" else None
            )

        if self.pairing == "paired" and r2 is None:
            raise ValueError("Paired sample requires R2 FASTQ")

        self.input_files = [r1] + ([r2] if r2 else [])

    @property
    def fastq_r1(self) -> Path:
        return self.input_files[0]

    @property
    def fastq_r2(self) -> Path | None:
        return self.input_files[1] if len(self.input_files) > 1 else None


@dataclass
class Genome:
    """Lightweight genome description"""

    species: str
    revision: int
    prebuild_prefix: str
    genetic_code: Any | None = None

    @property
    def base(self) -> Path:
        base_dir = Path("/machine/ffs/prebuild/externals/ppg2/clara/ensembl")
        return base_dir / self.species

    @property
    def name(self) -> str:
        return f"{self.species}_{self.revision}"

    @property
    def fasta(self) -> Path:
        return self.base / self.prebuild_prefix / "genome.fasta"

    @property
    def url(self) -> str:
        return f"{self.base}/{self.prebuild_prefix}/genome.url"

    def __repr__(self) -> str:
        return (
            f"Genome(species={self.species}, revision={self.revision}, "
            f"base={self.base}, fasta_file={self.fasta}, base_url={self.url})"
        )

    def __str__(self) -> str:
        return self.__repr__()


@dataclass
class HardFilterThresholds:
    """Thresholds for the hard-filter expression applied to tumour variants.

    Attributes
    ----------
    min_dp : int
        Minimum total read depth (``FMT/DP``).
    min_ad : int
        Minimum alt-allele depth (``FMT/AD[1]``).
    min_vaf : float
        Minimum variant allele frequency (``FMT/AD[1]/FMT/DP``).
    """

    min_dp: int = 10
    min_ad: int = 3
    min_vaf: float = 0.05

    def as_expression(self) -> str:
        """Return a bcftools ``-i`` filter expression for these thresholds.

        Returns
        -------
        str
            e.g. ``"FMT/DP>=10 && FMT/AD[1]>=3 && (FMT/AD[1]/FMT/DP)>=0.05"``
        """
        return (
            f"FMT/DP>={self.min_dp} && "
            f"FMT/AD[1]>={self.min_ad} && "
            f"(FMT/AD[1]/FMT/DP)>={self.min_vaf}"
        )
