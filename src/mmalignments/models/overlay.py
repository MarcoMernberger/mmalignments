from __future__ import annotations

from dataclasses import dataclass, field, fields, replace
from functools import cached_property
from pathlib import Path
from typing import (
    IO,
    Any,
    Generic,
    Iterator,
    Mapping,
    Self,  # # type: ignore[import] Python 3.11+
    TypeVar,
)

from mmalignments.models.resources import ResourceConfig
from mmalignments.models.tags import (
    Method,
    Omics,
    Stage,
    State,
)

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

    # def resolve(self, *patches: ElementTag | PartialElementTag | None) -> ElementTag:
    #     """
    #     Resolve the element tag by applying the given patches.

    #     Each patch can be a PartialElementTag, an ElementTag, or None. The
    #     patches are applied in order, with later patches overriding earlier
    #     ones.

    #     Returns
    #     -------
    #     ElementTag
    #         The resolved element tag.
    #     """
    #     result = self
    #     for patch in patches:
    #         result = result.patch(patch)
    #     return result

    def bump(
        self,
        root: str | None = None,
        level: int | None = None,
        stage: Stage | None = None,
        method: Method | None = None,
        state: State | None = None,
        omics: Omics | None = None,
        flag: str | None = None,
    ) -> ElementTag:
        """
        Create a new ElementTag by bumping the current tag with the provided
        values.

        Parameters
        ----------
        root : str | None, optional
            The new root value, by default None
        level : int | None, optional
            The new level value, by default None
        stage : Stage | None, optional
            The new stage value, by default None
        method : Method | None, optional
            The new method value, by default None
        state : State | None, optional
            The new state value, by default None
        omics : Omics | None, optional
            The new omics value, by default None
        flag : str | None, optional
            The new flag value, by default None

        Returns
        -------
        ElementTag
            A new ElementTag instance with the bumped values.
        """
        return ElementTag.from_prior(
            self,
            root=root or self.root,
            level=level if level is not None else self.level + 1,
            stage=stage or self.stage,
            method=method or self.method,
            state=state or self.state,
            omics=omics or self.omics,
            flag=flag or None,  # never inherit
        )

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
        prior = ElementTag(
            root=prior.root,
            level=prior.level,
            stage=prior.stage,
            method=prior.method,
            state=prior.state,
            omics=prior.omics,
            flag=None,
        )
        patched = prior.patch(PartialElementTag(flag=None)) # never inherit
        patched = patched.resolve(base, tag)
        # override with kwargs
        patched = patched.patch(PartialElementTag(**kwargs))
        return patched


################################################################################
# OutSpec
################################################################################


@dataclass(frozen=True)
class FileSpec:
    "Short specification for a single additional output file"

    stem: str
    ext: str


@dataclass(frozen=True)
class PartialOutSpec(Overlay["PartialOutSpec"]):
    """
    Partial specification for the output of an Element, allowing for optional
    overrides.
    """

    stem: str | None = None
    folder: Path | None = None
    prefixes: list[str] | None = None
    suffixes: list[str] | None = None
    ext: str | None = None
    exts: tuple[str, ...] | None = None
    additional_output: dict[str, FileSpec] | None = None


@dataclass(frozen=True)
class OutSpec(Overlay[PartialOutSpec]):
    """
    A specification for the output of an Element, including the output
    directory, file name(s), and extension(s). Used to generate default names
    and paths for output files, and to override defaults where needed.

    The primary output is the cross product of ``prefixes`` x ``suffixes``
    for ``ext`` (plus any ``exts``). In the common case (single prefix/suffix,
    both ``""``) this yields exactly one primary file. Setting e.g.
    ``suffixes=["_R1", "_R2"]`` yields two primary files sharing one stem
    (paired FASTQ). At most one of ``prefixes``/``suffixes`` may hold more
    than one entry at a time, to avoid accidental combinatorial explosions.

    ``stem`` may be left as ``None`` at construction time (e.g. when an
    ``OutSpec`` is declared statically, before an ``ElementTag`` exists) and
    resolved later via :meth:`resolve_stem`. At least one of ``stem``,
    ``prefixes``, or ``suffixes`` must contribute non-empty content, or no
    filename can be constructed.
    """

    stem: str | None = None
    folder: Path | None = None
    prefixes: list[str] = field(default_factory=lambda: [""])
    suffixes: list[str] = field(default_factory=lambda: [""])
    ext: str = "parquet"
    exts: tuple[str, ...] | None = None
    additional_output: dict[str, FileSpec] = field(
        default_factory=dict
    )  # files in addition to primary output, e.g. {"tsv": FileSpec(stem="results", ext="tsv")} # noqa: E501

    _REQUIRED = ("ext",)  # fields that must not be None (stem is allowed to be None)

    def __post_init__(self) -> None:
        missing = [f for f in self._REQUIRED if getattr(self, f) is None]
        if missing:
            raise ValueError(
                f"OutSpec: required field(s) must not be None: {', '.join(missing)}"  # noqa: E501
            )

        if not self.prefixes and not self.suffixes and not self.stem:
            raise ValueError(
                "OutSpec: cannot construct a filename - stem, prefixes, and suffixes are all empty. Set 'stem' (directly, or via resolve_stem()), or a non-empty prefix/suffix."  # noqa: E501
            )

        # if len(self.prefixes) > 1 and len(self.suffixes) > 1:
        #     raise ValueError(
        #         "OutSpec: only one of 'prefixes'/'suffixes' may have multiple "
        #         "entries at a time."
        #     )

    ############################################################################
    # Name / path construction
    ############################################################################

    def _make_path(self, filename: str) -> Path:
        if self.folder is None:
            return Path(filename)
        return self.folder / filename

    def _make_file_name(self, prefix: str, stem: str, suffix: str, ext: str) -> str:
        base = f"{prefix}{stem}{suffix}"
        if not base:
            raise ValueError(
                "OutSpec: cannot construct a filename - stem, prefix, and "
                "suffix are all empty. Set 'stem' (directly, or via "
                "resolve_stem()), or a non-empty prefix/suffix."
            )
        return f"{base}.{ext}"

    @cached_property
    def output_groups(self) -> dict[str, list[Path]]:
        """
        All output files, grouped by key. The primary extension (``self.ext``)
        and any ``exts`` each map to a list of paths (cross product of
        ``prefixes`` x ``suffixes``). Each ``additional_output`` entry maps to
        a single-element list.
        """
        stem = self.stem or ""
        groups: dict[str, list[Path]] = {}

        for ext in [self.ext, *(self.exts or ())]:
            group_paths = [
                self._make_path(
                    self._make_file_name(
                        prefix=prefix, stem=stem, suffix=suffix, ext=ext
                    )
                )
                for prefix in self.prefixes
                for suffix in self.suffixes
            ]
            groups[ext] = group_paths

        for key, file_spec in self.additional_output.items():
            fullname = self._make_file_name(
                prefix="", stem=file_spec.stem, suffix="", ext=file_spec.ext
            )
            groups[key] = [self._make_path(fullname)]

        return groups

    ############################################################################
    # Stem resolution
    ############################################################################

    def resolve_stem(self, default: str) -> OutSpec:
        """
        Return a new OutSpec with ``stem`` set to ``default`` if it is not
        already set. No-op if ``stem`` is already set.

        Parameters
        ----------
        default : str
            The stem to use if none is currently set.

        Returns
        -------
        OutSpec
            A new OutSpec instance with ``stem`` resolved.
        """
        if self.stem is not None:
            return self
        return replace(self, stem=default)

    ############################################################################
    # Mutation
    ############################################################################

    def add_output(self, key: str, file_spec: FileSpec) -> OutSpec:
        """
        Add an additional output file specification to the OutSpec.
        """
        new_additional_output = dict(self.additional_output)
        new_additional_output[key] = file_spec
        return replace(self, additional_output=new_additional_output)

    # def resolve(self, *patches: OutSpec | PartialOutSpec | None) -> OutSpec:
    #     """
    #     Resolve the OutSpec by applying the given patches, in order.
    #     """
    #     result = self
    #     for patch in patches:
    #         result = result.patch(patch)
    #     return result

    ############################################################################
    # Iteration / access
    ############################################################################

    def iter_grouped(self) -> Iterator[tuple[str, list[Path]]]:
        """
        Iterate over (key, paths) pairs, primary extension first.
        """
        yield self.ext, self.output_groups.get(self.ext, [])
        for key, paths in self.output_groups.items():
            if key == self.ext:
                continue
            yield key, paths

    def iterfiles(self) -> Iterator[tuple[str, Path]]:
        """
        Flatten iter_grouped() into (key, path) pairs, one per file.
        """
        for key, paths in self.iter_grouped():
            for path in paths:
                yield key, path

    @property
    def files(self) -> list[Path]:
        """
        Flat list of every path this OutSpec produces, across all groups
        (primary, exts, and additional_output combined).
        """
        return [path for paths in self.output_groups.values() for path in paths]

    def get(self, key: str) -> Path | tuple[Path, ...]:
        """
        Return the path(s) for a given key (a primary/exts extension, or an
        additional_output key). Returns a single Path if the group has
        exactly one file, otherwise a tuple of Paths (e.g. paired FASTQ).

        Raises
        ------
        KeyError
            If ``key`` is not a known output group.
        """
        try:
            paths = self.output_groups[key]
        except KeyError:
            raise KeyError(f"OutSpec has no output group '{key}'.") from None

        if len(paths) == 1:
            return paths[0]
        return tuple(paths)

    def map(self) -> dict[str, Path | tuple[Path, ...]]:
        """
        Mapping of key -> Path (single-file group) or tuple[Path, ...]
        (multi-file group, e.g. paired FASTQ). Covers every output group.
        """
        return {key: self.get(key) for key in self.output_groups}

    def filename(self) -> str:
        paths = self.output_groups.get(self.ext, [])
        if len(paths) == 0:
            raise ValueError("No output file defined for the primary extension.")
        if len(paths) > 1:
            raise ValueError(
                "Multiple output files defined for the primary extension; "
                "use 'output_groups' or 'iter_grouped()' instead of "
                "'filename()'/'path()'/'file'."
            )
        return paths[0].name

    def path(self, key: str | None = None) -> Path:
        if not key:
            filename = self.filename()
            return self._make_path(filename)
        else:
            paths = self.output_groups.get(key, [])
            if len(paths) == 0:
                raise ValueError(
                    f"No output file defined for the output group '{key}'."
                )  # noqa: E501
            if len(paths) > 1:
                raise ValueError(
                    f"Output group '{key}' has multiple files; use 'output_groups' or 'iter_grouped()' instead of 'path()'/'file'."  # noqa: E501
                )
            return paths[0]

    @property
    def file(self) -> Path:
        return self.path()

    ############################################################################
    # Construction helpers
    ############################################################################

    @classmethod
    def from_tag(
        cls, tag: ElementTag, *patches: OutSpec | PartialOutSpec | None, **kwargs
    ) -> OutSpec:
        """
        Create an OutSpec from an ElementTag, optionally applying a patch.
        """
        base = OutSpec(
            stem=tag.default_name,
            folder=Path("results") / tag.root,
            ext="tsv",
        )
        override = OutSpec(**kwargs) if kwargs else None
        return base.resolve(*patches, override)


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

    # def resolve(
    #     self, *patches: ExternalRunConfig | PartialExternalRunConfig | None
    # ) -> ExternalRunConfig:
    #     result = self
    #     for patch in patches:
    #         result = result.patch(patch)
    #     return result


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

    def patch(self, other: Overlay["Params"] | None):
        if other is None:
            return self

        return Params(
            **{
                **self._params,
                **other.to_dict(),
            }
        )

    def get(self, key, default=None):
        return self._params.get(key, default)

    # def resolve(self, *patches: Params | None):
    #     result = self

    #     for patch in patches:
    #         result = result.patch(patch)

    #     return result

    def override(self, **kwargs: Any) -> Params:
        """
        override convenience method to not break External.

        Returns
        -------
        Params
            A new Params instance with the overridden values.
        """
        return self.patch(Params(**kwargs))

    @cached_property
    def determinants(self) -> tuple[str, ...]:
        """
        Return a list of strings that uniquely identify the current Params instance.
        This can be used for caching or comparison purposes.

        Returns
        -------
        tuple[str]
            A tuple of strings representing the unique determinants of the Params.
        """
        return tuple(f"{k}={v}" for k, v in sorted(self._params.items()))


########################################################################################
# Shortcuts
########################################################################################

Tag = PartialElementTag
Out = PartialOutSpec
Par = Params
Ext = PartialExternalRunConfig

OutputSpec = OutSpec
Ecfg = ExternalRunConfig

TagType = ElementTag | PartialElementTag
OutType = OutSpec | PartialOutSpec
ParType = Params
CfgType = ExternalRunConfig | PartialExternalRunConfig


def from_prior(
    prior: ElementTag,
    tag: TagType | None = None,
    **kwargs,
) -> ElementTag:
    """
    Compatibility wrapper for modules that import from_prior from the overlay layer.
    """
    return ElementTag.from_prior(prior, tag=tag, **kwargs)


################################################################################
# StepValues
################################################################################

T = TypeVar("T")


class SpecValues(Generic[T]):
    """Per-step override container for grabfast-style multi-Element calls.
    Unset steps fall back to `default`."""

    def __init__(self, default: T | None, **kwargs) -> None:
        self.default = default
        for step in kwargs:
            setattr(self, step, kwargs.get(step, None))

    def get(self, step: str) -> T | None:
        return getattr(self, step, None) or self.default


def pick_spec(
    value: T | SpecValues[T] | None,
    step: str,
) -> T | None:
    if value is None:
        return None

    if isinstance(value, SpecValues):
        return value.get(step)

    return value
