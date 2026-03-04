from dataclasses import dataclass
from enum import Enum
from functools import cached_property


class Omics(str, Enum):
    DNA = "dna"
    RNA = "rna"
    PROTEIN = "protein"
    METABOLOME = "metabolome"


class Stage(str, Enum):
    QC = "qc"
    PREP = "prep"
    ALIGN = "align"
    CALL = "call"
    INTEGRATE = "integrate"
    QUANT = "quant"
    DIFF = "diff"


class State(str, Enum):
    RAW = "raw"
    TRIM = "trim"
    SORT = "sort"
    FILTER = "filter"
    NORMAL = "normal"
    DEDUP = "dedup"
    RECAL = "recal"
    HARMONIZED = "harmonized"
    INDEX = "indexed"


class Method(str, Enum):
    FASTP = "fastp"
    BWAMEM2 = "bwamem2"
    GATK = "gatk"
    STAR = "star"
    SALMON = "salmon"
    DESEQ2 = "deseq2"
    FASTQC = "fastqc"
    MOSDEPTH = "mosdepth"
    CONTIGMAP = "contigmap"
    BEDTOOLS = "bedtools"
    SAMTOOLS = "samtools"
    PICARD = "picard"


def level(level: int) -> str:
    return f"{level:02d}"


@dataclass(frozen=True)
class ElementTag:
    root: str  # sample, e.g. Kidney_1
    level: int  # No. in pipeline chain
    stage: Stage  # e.g. align
    method: Method  # e.g. bwamem2
    state: State  # e.g. raw
    omics: Omics | None  # Omics domain
    ext: str | None = None  # file extension, e.g. tsv
    param: str | None = None  # parameter string, e.g. "minimap2-preset=sr"

    @cached_property
    def default_name(self) -> str:
        to_join = [
            self.root,
            level(self.level),
        ]
        if self.omics:
            to_join.append(
                self.omics,
            )
        to_join.extend(
            [
                self.stage,
                self.method,
                self.state,
            ]
        )
        if self.param:
            to_join.append(self.param)
        return ".".join(to_join)

    @cached_property
    def default_output(self) -> str:
        return self.default_name + (f".{self.ext}" if self.ext else "")
