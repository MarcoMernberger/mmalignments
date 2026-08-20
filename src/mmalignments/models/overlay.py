from __future__ import annotations

from dataclasses import dataclass, field, fields, replace
from functools import cached_property
from pathlib import Path
from typing import Any, IO, Mapping, Iterator, Generic, TypeVar
from mmalignments.models.resources import ResourceConfig
from mmalignments.models.tags import (
    Method,
    Omics,
    Stage,
    State,
)

from typing import Self  # # type: ignore[import] Python 3.11+

TPartial = TypeVar("TPartial", bound="Overlay")


@dataclass(frozen=True)
class Overlay(Generic[TPartial]):
    """
    Base class for immutable configuration overlays.
    """

    def to_dict(self, drop_none: bool = True) -> dict[str, Any]:
        """
        Convert the overlay to a dictionary.

        This method iterates over all fields of the dataclass and constructs a
        dictionary representation. If `drop_none` is True, fields with a value of None
        are omitted.

        Parameters
        ----------
        drop_none : bool, optional
            Whether to drop fields with a value of None, by default True

        Returns
        -------
        dict[str, Any]
            A dictionary representation of the overlay.
        """
        result = {}

        for f in fields(self):
            value = getattr(self, f.name)

            if drop_none and value is None:
                continue

            result[f.name] = value

        return result

    def items(self):
        """Return an iterable of (field_name, field_value) pairs."""
        return self.to_dict(drop_none=False).items()

    def keys(self):
        """Return an iterable of field names."""
        return self.to_dict(drop_none=False).keys()

    def values(self):
        """Return an iterable of field values."""
        return self.to_dict(drop_none=False).values()

    def __getitem__(self, key: str):
        """Allow dictionary-like access to fields."""
        return getattr(self, key)

    def __iter__(self):
        """Return an iterator over the field names."""
        return iter(self.keys())

    def __contains__(self, key: str) -> bool:
        """Check if a field name exists in the overlay."""
        return key in self.keys()

    def update(self, **kwargs):
        """
        Override fields of the overlay with the provided keyword arguments using
        dataclass replace.

        Returns
        -------
        Overlay
            A new overlay instance with the overridden fields.
        """
        return replace(self, **kwargs)

    def patch(self, other: Self | TPartial | None) -> Self:
        """
        Patch the current overlay with another overlay.

        This method creates a new overlay instance where fields from the `other`
        overlay override the corresponding fields in the current overlay. If `other`
        is None, the current overlay is returned unchanged.

        Parameters
        ----------
        other : Self | TPartial | None
            The overlay instance to patch with, or None.

        Returns
        -------
        Overlay
            A new overlay instance with the patched fields.
        """
        if other is None:
            return self

        values = other.to_dict()

        return replace(self, **values)

    # def resolve(self, *patches: Overlay | None) -> Overlay:
    #     raise NotImplementedError(
    #         "Subclasses must implement the resolve method to apply patches."
    #     )

    def resolve(self, *patches: Self | TPartial | None) -> Self:
        result = self
        for patch in patches:
            result = result.patch(patch)
        return result


################################################################################
# Tags
################################################################################


@dataclass(frozen=True)
class PartialElementTag(Overlay):
    """Partial ElementTag, allowing for optional overrides of the ElementTag fields."""

    root: str | None = None
    level: int | None = None
    stage: Stage | None = None
    method: Method | None = None
    state: State | None = None
    omics: Omics | None = None
    flag: str | None = None


@dataclass(frozen=True)
class ElementTag(Overlay[PartialElementTag]):
    """
    ElementTag tags element in the overlay, describing the current step in the pipeline.

    This class encapsulates the various attributes that define an element tag,
    including its root, level, stage, method, state, omics, extension, and parameters.
    """

    root: str  # a short name for the sample, e.g. "sample1"
    level: int  # the level of the element in the execution order of the pipeline
    stage: Stage  # the stage of the pipeline, e.g. "preprocess", "align", "postprocess"
    method: (
        Method  # the method used for the current step, e.g. "bwa", "bowtie2", "star"
    )
    state: State  # the state of the element, e.g. "raw", "trimmed", "aligned", "sorted"
    omics: Omics | None = None  # the omics type, e.g. "DNA", "RNA", "proteome"
    flag: str | None = (
        None  # optional parameter flags, e.g. "trimmed", "deduped", "filtered"
    )

    _REQUIRED = (
        "root",
        "level",
        "stage",
        "method",
        "state",
    )  # fields that must not be None

    def __post_init__(self) -> None:
        """
        Validate that required fields are not None.

        Raises
        ------
        ValueError
            If any of the required fields are None.
        """
        missing = [f for f in self._REQUIRED if getattr(self, f) is None]
        if missing:
            raise ValueError(
                f"ElementTag: required field(s) must not be None: {', '.join(missing)}"  # noqa: E501
            )

    def format_level(self, level: int) -> str:
        return f"S{level:02d}"

    @cached_property
    def default_name(self) -> str:
        """
        Get the default name for the element tag, to be used for default output file
        names.
        The default name is constructed by joining the root, level, omics (if present),
        stage, method, state, and param (if present) with dots.

        Returns
        -------
        str
            The default name for the element tag.
        """
        to_join = [
            self.root,
            self.format_level(self.level),
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
        if self.flag:
            to_join.append(self.flag)
        return ".".join(to_join)

    def resolve(self, *patches: PartialElementTag | ElementTag | None) -> ElementTag:
        """
        Resolve the element tag by applying the given patches.

        Each patch can be a PartialElementTag, an ElementTag, or None. The patches
        are applied in order, with later patches overriding earlier ones.

        Returns
        -------
        ElementTag
            The resolved element tag.
        """
        result = self
        for patch in patches:
            result = result.patch(patch)
        return result

    @classmethod
    def from_prior(
        cls,
        prior: ElementTag,
        tag: PartialElementTag | ElementTag | None = None,
        **kwargs,
    ) -> ElementTag:
        """Return a new ElementTag by deriving fields from a prior ElementTag, from
        the previous Element in chain. Changed fields for the new ElementTag can be
        specified as keyword arguments.

        This lets callers specify only the fields that differ from the computed
        default::

            ensemblmap(element, tag=PartialElementTag(root="sample"))
        """
        base = PartialElementTag(
            root=prior.root,  # default to same root as prior since usually same sample
            level=prior.level
            + 1,  # default to one level higher than prior since usually is  next step
            stage=prior.stage,  # default to same stage as prior
            state=None,  # default to None since state often changes at each step
            method=prior.method,  # default to same method as prior
            omics=prior.omics,  # default to same omics as prior
            flag=None,  # flag are unlikely to be inherited, so default to None
        )

        patched = prior.resolve(base, tag)
        # override with kwargs
        patched = patched.patch(PartialElementTag(**kwargs))
        return patched


################################################################################
# OutputSpec
################################################################################


@dataclass(frozen=True)
class FileSpec:
    "Short specification for an additional output file"

    name: str
    ext: str


@dataclass(frozen=True)
class PartialOutputSpec(Overlay["PartialOutputSpec"]):
    """
    Partial specification for the output of an Element, allowing for optional overrides.
    """

    stem: str | None = None
    outdir: Path | None = None
    prefix: str | None = None
    suffix: str | None = None
    ext: str | None = None
    exts: tuple[str, ...] | None = None
    additional_output: dict[str, FileSpec] | None = None


@dataclass(frozen=True)
class OutputSpec(Overlay[PartialOutputSpec]):
    """
    A specification for the output of an Element, including the output directory,
    file name, and extension. Used to generate default names and paths for output files,
    and to override defaults specifically.
    """

    stem: str | None = None
    outdir: Path | None = None
    prefix: str = ""
    suffix: str = ""
    ext: str = "parquet"
    exts: tuple[str, ...] | None = None
    additional_output: dict[str, FileSpec] = field(
        default_factory=dict
    )  # files in addition to primary output, e.g. {"tsv": FileSpec(name="results", ext="tsv")} # noqa: E501

    _REQUIRED = (
        "stem",
        "ext",
        "prefix",
        "suffix",
    )  # fields that must not be None

    def __post_init__(self) -> None:
        """
        Validate that required fields are not None.

        Raises
        ------
        ValueError
            If any of the required fields are None.
        """
        missing = [f for f in self._REQUIRED if getattr(self, f) is None]
        if missing:
            raise ValueError(
                f"OutputSpec: required field(s) must not be None: {', '.join(missing)}"  # noqa: E501
            )

    def add_output(self, key: str, file_spec: FileSpec) -> OutputSpec:
        """
        Add an additional output file specification to the OutputSpec.

        Parameters
        ----------
        key : str
            The key for the additional output file.
        file_spec : FileSpec
            The specification for the additional output file.

        Returns
        -------
        OutputSpec
            A new OutputSpec instance with the additional output file added.
        """
        new_additional_output = dict(self.additional_output)
        new_additional_output[key] = file_spec
        return replace(self, additional_output=new_additional_output)

    def resolve(self, *patches: OutputSpec | PartialOutputSpec | None) -> OutputSpec:
        """
        Resolve the OutputSpec by applying the given patches.

        Each patch can be an OutputSpec, PartialOutputSpec, or None. The patches are
        applied in order, with later patches overriding earlier ones.

        Returns
        -------
        OutputSpec
            The resolved OutputSpec.
        """
        result = self
        for patch in patches:
            result = result.patch(patch)
        return result

    @property
    def files(self) -> dict[str, Path]:
        return dict(self.iterfiles())

    @property
    def file(self) -> Path:
        return self.path()

    def iterfiles(self, default_name: str | None = None) -> Iterator[tuple[str, Path]]:
        name = self.stem or default_name
        if name is None:
            return iter(())

        # primary
        yield (
            self.ext,
            self._make_path(self._make_file_name(name, self.ext)),
        )
        # primary additional extensions
        for ext in self.exts or ():
            yield (
                ext,
                self._make_path(self._make_file_name(name, ext)),
            )
        # additional named files
        for key, extra in self.additional_output.items():
            yield (
                key,
                self._make_path(self._make_file_name(extra.name, extra.ext)),
            )

    def _make_path(self, filename: str) -> Path:
        if self.outdir is None:
            return Path(filename)

        return self.outdir / filename

    def _make_file_name(self, name: str, ext: str) -> str:
        return f"{self.prefix}{name}{self.suffix}.{ext}"

    def filename(self, default_stem: str | None = None) -> str:
        name = self.stem or default_stem
        if name is None:
            raise ValueError("No name provided for output file.")
        return self._make_file_name(name, self.ext)

    def path(self, default_stem: str | None = None) -> Path:
        filename = self.filename(default_stem=default_stem)
        return self._make_path(filename)


########################################################################################
# ExternalRunConfig
########################################################################################


@dataclass(frozen=True)
class PartialExternalRunConfig(Overlay["PartialExternalRunConfig"]):
    """Partial configuration for running an External tool."""

    cwd: Path | None = None
    env: dict[str, str] | None = None
    pipe_output: bool | None = None
    check: bool | None = None
    timeout: float | None = None
    threads: int | None = None
    multi: bool | None = None
    stdout: Path | None | IO = None
    stderr: Path | None | IO = None
    append: bool | None = None
    log_dir: Path | None = None


@dataclass(frozen=True)
class ExternalRunConfig(Overlay[PartialExternalRunConfig]):
    """Configuration for running an External tool."""

    cwd: Path | None = None
    env: dict[str, str] | None = None
    pipe_output: bool = True
    check: bool = True
    timeout: float | None = None
    threads: int = field(
        default_factory=lambda: ResourceConfig.detect().threads
    )  # noqa: E501
    multi: bool = True
    stdout: Path | None | IO = None
    stderr: Path | None | IO = None
    append: bool = False
    log_dir: Path | None = None

    def __post_init__(self) -> None:
        if self.stdout and self.stderr and self.stdout == self.stderr:
            raise ValueError(
                "stdout and stderr cannot be the same path or file object"
            )  # noqa: E501
        if self.log_dir and not isinstance(self.log_dir, Path):
            raise ValueError(f"log_dir must be a directory: {type(self.log_dir)}")
        if self.threads < 1:
            raise ValueError("threads must be a positive integer")

    def resolve(
        self, *patches: ExternalRunConfig | PartialExternalRunConfig | None
    ) -> ExternalRunConfig:
        result = self
        for patch in patches:
            result = result.patch(patch)
        return result


################################################################################
# Params
################################################################################


@dataclass(frozen=True)
class Params(Overlay["Params"], Mapping[str, Any]):
    _params: dict[str, Any] = field(default_factory=dict)

    def __init__(self, **kwargs: Any):
        object.__setattr__(
            self,
            "_params",
            dict(kwargs),
        )

    def __getitem__(self, key: str):
        return self._params[key]

    def __iter__(self):
        return iter(self._params)

    def __len__(self):
        return len(self._params)

    def to_dict(self, drop_none=True):
        if drop_none:
            return {k: v for k, v in self._params.items() if v is not None}

        return dict(self._params)

    def patch(self, other: Params | None):
        if other is None:
            return self

        return Params(
            **{
                **self._params,
                **other._params,
            }
        )

    def get(self, key, default=None):
        return self._params.get(key, default)

    def resolve(self, *patches: Params | None):
        result = self

        for patch in patches:
            result = result.patch(patch)

        return result


########################################################################################
# Shortcuts
########################################################################################

Tag = PartialElementTag
Out = PartialOutputSpec
Par = Params
Ext = PartialExternalRunConfig

TagType = ElementTag | PartialElementTag
OutType = OutputSpec | PartialOutputSpec
ParType = Params
ExtType = ExternalRunConfig | PartialExternalRunConfig


def from_prior(
    prior: ElementTag,
    tag: TagType | None = None,
    **kwargs,
) -> ElementTag:
    """
    Compatibility wrapper for modules that import from_prior from the overlay layer.
    """
    return ElementTag.from_prior(prior, tag=tag, **kwargs)
