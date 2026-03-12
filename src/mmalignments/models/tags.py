import dataclasses
from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from deprecated import deprecated  # type: ignore
from typing import Any


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
    INPUT = "input"


class State(str, Enum):
    RAW = "raw"
    TRIM = "trim"
    SORT = "sort"
    MAP = "map"
    FILTER = "filter"
    NORMAL = "normal"
    COUNT = "count"
    DEDUP = "dedup"
    RECAL = "recal"
    MODEL = "model"
    HARMONIZED = "harmonic"
    INDEX = "index"
    PADDED = "pad"
    ANNOTATE = "annotate"
    MERGED = "merge"
    PILE = "pile"
    STAT = "stat"
    LOCI = "loci"
    METRIC = "metric"
    REPORT = "report"


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
    MULTIQC = "multiqc"
    CHECK = "check"
    BCFTOOLS = "bcftools"
    CUSTOM = "custom"  # for custom python functions


def level(level: int) -> str:
    return f"S{level:02d}"


@dataclass(frozen=True)
class ElementTag:
    root: str  # sample, e.g. Kidney_1
    level: int  # No. in pipeline chain
    stage: Stage  # e.g. align
    method: Method  # e.g. bwamem2
    state: State  # e.g. raw
    omics: Omics | None = None  # Omics domain
    ext: str | None = None  # file extension, e.g. tsv
    param: str | None = None  # parameter string, e.g. "minimap2-preset=sr"

    _REQUIRED = (
        "root",
        "level",
        "stage",
        "method",
        "state",
    )  # fields that must not be None

    def __post_init__(self) -> None:
        missing = [f for f in self._REQUIRED if getattr(self, f) is None]
        if missing:
            raise ValueError(
                f"ElementTag: required field(s) must not be None: {', '.join(missing)}"
            )

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

    def __getitem__(self, key: str) -> Any:
        """Allows tag["state"] access to the fields."""
        if key in [f.name for f in dataclasses.fields(self)]:
            return getattr(self, key)
        raise KeyError(f"{key} is not a valid PartialElementTag field")

    def keys(self):
        return [f.name for f in dataclasses.fields(self)]

    def values(self):
        return [getattr(self, f.name) for f in dataclasses.fields(self)]

    def items(self):
        return [(f.name, getattr(self, f.name)) for f in dataclasses.fields(self)]

    def merge(self, override: "PartialElementTag | ElementTag | None") -> "ElementTag":
        """
        Return a new ElementTag using *self* as default,
        overriding with any non-None field from *override*.

        This lets callers specify only the fields that differ from the computed
        default.


        Parameters
        ----------
        patch : PartialElementTag | None
            A PartialElementTag containing fields to override in the base.

        Returns
        -------
        ElementTag
            A fully initialized ElementTag instance.

        Examples
        --------

            default_tag.merge(tag)
        """
        return PartialElementTag(**self).merge(override).resolve()


@dataclass(frozen=True, slots=True)
class PartialElementTag:
    """Partial ElementTag for use as a patch argument to merge_tag.

    All fields are optional so callers can specify only the fields that should
    differ from the computed default::

        ensemblmap(element, tag=PartialElementTag(root="sample"))
    """

    root: str | None = None
    level: int | None = None
    stage: Stage | None = None
    method: Method | None = None
    state: State | None = None
    omics: Omics | None = None
    ext: str | None = None
    param: str | None = None

    def resolve(self) -> ElementTag:
        """
        resolve Creates a full ElementTag from this PartialElementTag.
        ElementTag.__post_init__ ensures that all required fields are set.

        Returns
        -------
        ElementTag
            A fully initialized ElementTag instance.
        """
        return ElementTag(**dataclasses.asdict(self))

    def patch(self, **overrides) -> "PartialElementTag":
        """
        patch patches this PartialElementTag with any non-None fields from
        overrides, returning a new PartialElementTag.

        Returns
        -------
        PartialElementTag
            A new PartialElementTag instance with the applied overrides.
        """
        data = dataclasses.asdict(self)
        data.update({k: v for k, v in overrides.items() if v is not None})
        return PartialElementTag(**data)

    def merge(
        self, other: "PartialElementTag | ElementTag | None"
    ) -> "PartialElementTag":
        """
        Return a new PartialElementTag using self as default,
        overriding with any non-None field from *other*.

        This lets callers specify only the fields that differ from the computed
        default.

        Parameters
        ----------
        other : PartialElementTag | None
            A PartialElementTag containing fields to override in the base.

        Returns
        -------
        PartialElementTag
            A new PartialElementTag instance with the applied overrides.
        """
        if other is None:
            return self
        return self.patch(
            **{k: v for k, v in dataclasses.asdict(other).items() if v is not None}
        )

    def __getitem__(self, key: str) -> Any:
        """Allows tag["state"] access to the fields."""
        if key in [f.name for f in dataclasses.fields(self)]:
            return getattr(self, key)
        raise KeyError(f"{key} is not a valid PartialElementTag field")

    def keys(self):
        return [f.name for f in dataclasses.fields(self)]

    def values(self):
        return [getattr(self, f.name) for f in dataclasses.fields(self)]

    def items(self):
        return [(f.name, getattr(self, f.name)) for f in dataclasses.fields(self)]


def Tag(tag: PartialElementTag | None = None, **kwargs) -> PartialElementTag:
    """
    Convenience function to create a PartialElementTag, optionally patching an
    existing tag.

    Returns
    -------
    PartialElementTag
        A partially initialized ElementTag instance.
    """
    base = tag or PartialElementTag()
    return base.patch(**kwargs)


def from_prior(
    prior: ElementTag, tag: PartialElementTag | ElementTag | None = None, **kwargs
) -> ElementTag:
    """Return a new ElementTag by deriving fields from a prior ElementTag, from
    the previous Element in chain. Changed fields for the new ElementTag can be
    specified as keyword arguments.

    This lets callers specify only the fields that differ from the computed
    default::

        ensemblmap(element, tag=PartialElementTag(root="sample"))
    """
    base = PartialElementTag(
        root=prior.root,  #   default to same root as prior since usually the same sample, but can be overridden if needed
        level=prior.level
        + 1,  # default to one level higher than prior since usually the next step in the pipeline, but can be overridden if needed
        stage=prior.stage,  # default to same stage as prior since often the same stage, but can be overridden if needed
        state=None,  # default to None since state often changes at each step
        method=prior.method,  # default to same method as prior since often the same method, but can be overridden if needed
        omics=prior.omics,  # default to same omics as prior since usually the same sample, but can be overridden if needed
        ext=prior.ext,  # extension may or may not be inherited, but default to same as prior if not specified
        param=None,  # parameter flags are unlikely to be inherited, so default to None
    )
    # override with kwargs
    patched = base.patch(**kwargs)
    # override again with user-supplied tag
    if tag:
        patched = patched.merge(tag)

    return patched.resolve()


@deprecated(reason="Use from_prior(input_element.tag, **tag) instead")
def merge_tag(base: ElementTag, override: PartialElementTag | None) -> ElementTag:
    """
    Return a new ElementTag using *base* as default,
    overriding with any non-None field from *override*.

    This lets callers specify only the fields that differ from the computed
    default.


    Parameters
    ----------
    base : ElementTag | PartialElementTag
        The base ElementTag or PartialElementTag to use as default.
    patch : PartialElementTag | None
        A PartialElementTag containing fields to override in the base.

    Returns
    -------
    ElementTag
        A fully initialized ElementTag instance.

    Examples
    --------

        ensemblmap(element, tag=PartialElementTag(root="sample"))
    """
    return PartialElementTag(**base).merge(override).resolve()
