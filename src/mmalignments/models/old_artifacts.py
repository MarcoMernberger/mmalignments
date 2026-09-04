# from __future__ import annotations

# from abc import ABC
# from dataclasses import dataclass, field, fields, replace
# from functools import cached_property
# from pathlib import Path
# from types import MappingProxyType
# from typing import Any, Iterator, Literal, Mapping

# from pandas import DataFrame  # type: ignore[import]

# from mmalignments.core.annotations import ColumnSchema, View
# from mmalignments.models.tags import ElementTag
# from mmalignments.services.dependencies import (
#     DynamicValue,
#     file_signature,
#     stable_hash,
# )
# from mmalignments.services.io import read_frame, read_schema

# ###############################################################################
# # External Wrapper
# ###############################################################################

# ArtifactType = Literal["file", "value"]
# ArtifactLifeTime = Literal["persistent", "transient", "ephemeral"]


# @dataclass(frozen=True)
# class FileSpec:
#     name: str
#     ext: str


# @dataclass(frozen=True)
# class OutputSpec:
#     stem: str | None = None
#     outdir: Path | None = None
#     prefix: str = ""
#     suffix: str = ""
#     ext: str = "parquet"
#     exts: tuple[str, ...] | None = None
#     additional_output: dict[str, FileSpec] = field(
#         default_factory=dict
#     )  # files in addition to primary output, e.g. {"tsv": FileSpec(name="results", ext="tsv")}

#     def __post_init__(self):
#         extensions = set((self.ext,) + (self.exts or ()))
#         keys = set(self.additional_output.keys())
#         duplicate = keys & extensions
#         if duplicate:
#             raise ValueError(
#                 f"Duplicate keys in additional_output, conflicting with extensions: {duplicate}"
#             )

#     def _make_path(self, filename: str) -> Path:
#         if self.outdir is None:
#             return Path(filename)

#         return self.outdir / filename

#     def _make_file_name(self, name: str, ext: str) -> str:
#         return f"{self.prefix}{name}{self.suffix}.{ext}"

#     @property
#     def files(self) -> dict[str, Path]:
#         return dict(self.iterfiles())

#     def iterfiles(self, default_name: str | None = None) -> Iterator[tuple[str, Path]]:
#         name = self.stem or default_name
#         if name is None:
#             return

#         # primary
#         yield (
#             self.ext,
#             self._make_path(self._make_file_name(name, self.ext)),
#         )
#         # primary additional extensions
#         for ext in self.exts or ():
#             yield (
#                 ext,
#                 self._make_path(self._make_file_name(name, ext)),
#             )
#         # additional named files
#         for key, extra in self.additional_output.items():
#             yield (
#                 key,
#                 self._make_path(self._make_file_name(extra.name, extra.ext)),
#             )

#     def filename(self, default_name: str | None = None) -> str:
#         name = self.stem or default_name
#         if name is None:
#             raise ValueError("No name provided for output file.")
#         return self._make_file_name(name, self.ext)

#     def path(self, default_name: str | None = None) -> Path:
#         filename = self.filename(default_name=default_name)
#         if filename is None:
#             raise ValueError("No name provided for output file.")
#         return self._make_path(filename)

#     def merge(self, other: "OutputSpec | None") -> "OutputSpec":
#         if other is None:
#             return self

#         values = {}

#         for f in fields(self):
#             old = getattr(self, f.name)
#             new = getattr(other, f.name)

#             # keep old value if new is None
#             if new is None:
#                 values[f.name] = old

#             # merge dictionaries
#             elif isinstance(old, dict) and isinstance(new, dict):
#                 values[f.name] = {**old, **new}

#             # overwrite normal values
#             else:
#                 values[f.name] = new

#         return replace(self, **values)

#     def add_output(self, key: str, spec: FileSpec) -> OutputSpec:
#         if key == self.ext or key in (self.exts or ()):
#             raise ValueError(f"Output key '{key}' conflicts with extension '{key}'.")
#         return replace(
#             self,
#             additional_output={
#                 **self.additional_output,
#                 key: spec,
#             },
#         )


# @dataclass(frozen=True)
# class ArtifactDef:
#     cls: type[Artifact]
#     kwargs: MappingProxyType[str, Any] = field(
#         default_factory=lambda: MappingProxyType({})
#     )


# class ArtifactSet(Mapping):
#     __slots__ = ("_primary", "_primary_name", "_extras")

#     def __init__(
#         self, primary: Any, /, primary_name: str | None = None, **extras: Any
#     ):  # noqa: E501
#         if "primary" in extras:
#             raise KeyError("'primary' is reserved for the primary artifact.")
#         if primary_name is not None and primary_name in extras:
#             raise KeyError(
#                 f"'{primary_name}' collides with an extra artifact key."
#             )  # noqa: E501
#         object.__setattr__(self, "_primary", primary)
#         object.__setattr__(self, "_primary_name", primary_name)
#         object.__setattr__(self, "_extras", MappingProxyType(dict(extras)))

#     def __getitem__(self, key: str) -> Any:
#         if key == "primary" or (
#             self._primary_name is not None and key == self._primary_name
#         ):
#             return self._primary
#         return self._extras[key]

#     def __contains__(self, key: object) -> bool:
#         return (
#             key == "primary" or key == self._primary_name or key in self._extras
#         )  # noqa: E501

#     def __iter__(self) -> Iterator[str]:
#         yield (
#             self._primary_name if self._primary_name is not None else "primary"
#         )  # noqa: E501
#         yield from self._extras

#     def iterfiles(self) -> Iterator[Path]:
#         for a in self.values():
#             if hasattr(a, "resolve"):
#                 yield a.resolve()
#             elif isinstance(a, FastqArtifact):
#                 yield from a

#     def __len__(self) -> int:
#         return 1 + len(self._extras)

#     @property
#     def primary(self) -> Any:
#         return self._primary

#     @property
#     def primary_name(self) -> str | None:
#         return self._primary_name

#     @property
#     def outputs(self) -> tuple[Path, ...]:
#         return tuple(self.iterfiles())

#     def with_extra(self, key: str, value: Any) -> ArtifactSet:
#         if key == "primary" or key == self._primary_name or key in self._extras:
#             raise KeyError(f"Artifact '{key}' already exists.")
#         return ArtifactSet(
#             self._primary,
#             primary_name=self._primary_name,
#             **self._extras,
#             **{key: value},
#         )

#     def __repr__(self) -> str:
#         extras = ", ".join(f"{k}={v!r}" for k, v in self._extras.items())
#         name = self._primary_name or "primary"
#         return f"ArtifactSet({name}={self._primary!r}{', ' + extras if extras else ''})"  # noqa: E501

#     @classmethod
#     def from_any(
#         cls, x: ArtifactSet | Mapping[str, Any], primary_key: str = "primary"
#     ) -> ArtifactSet:
#         """
#         Convert a mapping to an ArtifactSet for backwards compatibility. If
#         the mapping does not contain a key for the primary artifact, a fallback
#         is attempted using the first available value or a known extension.
#         """
#         if isinstance(x, ArtifactSet):
#             return x
#         if isinstance(x, Mapping):
#             if primary_key not in x:
#                 if "parquet" in x:
#                     primary_name = "parquet"
#                     primary = x[primary_name]
#                 elif "tsv" in x:
#                     primary_name = "tsv"
#                     primary = x[primary_name]
#                 else:
#                     primary_name, primary = next(iter(x.items()), (None, None))
#                 if primary is None:
#                     raise KeyError(
#                         f"Expected key '{primary_key}' in mapping. No fallback found for primary artifact in {list(x.keys())}."  # noqa: E501
#                     )
#             else:
#                 primary = x[primary_key]
#                 primary_name = primary_key

#             return cls(
#                 primary,
#                 primary_name=primary_name,
#                 **{k: v for k, v in x.items() if k != primary_name},
#             )
#         raise TypeError(f"Cannot build ArtifactSet from {type(x)!r}")

#     @classmethod
#     def generate_file_artifacts(
#         cls,
#         tag: ElementTag,
#         spec: OutputSpec | None = None,
#         default_dir: Path | None = None,
#         index_column: str | None = None,
#     ) -> tuple[ArtifactSet, list[Path]]:

#         definitions = {
#             "primary": ArtifactDef(TableArtifact, {"index_column": index_column})
#         }
#         artifacts = cls.generate(tag, default_dir, spec, definitions=definitions)
#         outputs = list(artifacts.output_files().values())
#         return artifacts, outputs  # backwards compatibility: remove soon

#     @classmethod
#     def generate(
#         cls,
#         tag: ElementTag,
#         default_dir: Path | None = None,
#         spec: OutputSpec | None = None,
#         definitions: Mapping[str, ArtifactDef] | None = None,
#     ) -> ArtifactSet:
#         # ensure we have a name and directory
#         default = ArtifactDef(
#             TableArtifact
#         )  # assuming we have mostly tables, should be paths
#         definitions = definitions or {}
#         default_dir = default_dir or Path("results")
#         spec = OutputSpec(stem=tag.default_name, outdir=default_dir).merge(spec)
#         additional_artifacts = {}
#         primary_name = spec.ext or "primary"
#         primary = None
#         for key, path in spec.iterfiles(default_name=tag.default_name):
#             definition = definitions.get(key, default)
#             if key == primary_name:
#                 primary = definition.cls(
#                     path,
#                     **definition.kwargs,
#                 )
#             else:
#                 additional_artifacts[key] = definition.cls(
#                     path,
#                     **definition.kwargs,
#                 )
#         return ArtifactSet(
#             primary,
#             primary_name=primary_name,
#             **additional_artifacts,
#         )

#     def with_extras(self, extras: Mapping[str, Any]) -> "ArtifactSet":
#         """
#         Return a new ArtifactSet with additional artifacts.

#         Raises
#         ------
#         KeyError
#             If any artifact key already exists.
#         """
#         conflicts = set(extras) & set(self)

#         if conflicts:
#             raise KeyError(f"Artifact(s) already exist: {sorted(conflicts)}")

#         return ArtifactSet(
#             self._primary,
#             primary_name=self._primary_name,
#             **self._extras,
#             **dict(extras),
#         )

#     def merge(self, other: "ArtifactSet") -> "ArtifactSet":
#         """
#         Return a new ArtifactSet containing artifacts from both sets.

#         The primary artifact of `self` is preserved.
#         The primary artifact of `other` is added as an extra artifact if
#         it has a different name.

#         Raises
#         ------
#         KeyError
#             If artifact names collide.
#         """
#         if not isinstance(other, ArtifactSet):
#             other = ArtifactSet.from_any(other)

#         merged = dict(self._extras)

#         # add other's primary
#         other_primary_name = (
#             other._primary_name if other._primary_name is not None else "primary"
#         )

#         if other_primary_name != (self._primary_name or "primary"):
#             if other_primary_name in self:
#                 raise KeyError(f"Artifact '{other_primary_name}' already exists.")
#             merged[other_primary_name] = other._primary

#         # add other's extras
#         conflicts = set(merged) & set(other._extras)
#         if conflicts:
#             raise KeyError(f"Artifact(s) already exist: {sorted(conflicts)}")

#         merged.update(other._extras)

#         return ArtifactSet(
#             self._primary,
#             primary_name=self._primary_name,
#             **merged,
#         )

#     def output_files(self) -> Mapping[str, Path]:
#         return {
#             k: a.resolve() for k, a in self.items() if hasattr(a, "resolve")
#         }  # noqa: E501

#     def files(self) -> Iterator[Path]:
#         for a in self.values():
#             if hasattr(a, "resolve"):
#                 yield a.resolve()


# class Artifact(ABC):

#     def signature(self) -> str:
#         raise NotImplementedError


# @dataclass(frozen=True)
# class FastqArtifact(Artifact):
#     r1: Path
#     r2: Path | None = None
#     FASTQ_SUFFIXES: tuple[str, ...] = ("fq", "fastq", "fq.gz", "fastq.gz")

#     @property
#     def stem(self) -> str:
#         stem = self.r1.name
#         for suffix in self.FASTQ_SUFFIXES:
#             rem = "." + suffix
#             if stem.endswith(rem):
#                 stem = stem.removesuffix(rem)
#                 break
#         return stem

#     def signature(self) -> str:
#         for_hash = (
#             (file_signature(self.r1), file_signature(self.r2))
#             if self.r2
#             else (file_signature(self.r1),)
#         )
#         return stable_hash(for_hash)

#     @property
#     def paired(self) -> bool:
#         return self.r2 is not None

#     @property
#     def files(self) -> tuple[Path, ...]:
#         if self.paired:
#             return (self.r1.resolve(), self.r2.resolve())  # type: ignore[return-value]
#         else:
#             return (self.r1.resolve(),)

#     def __iter__(self) -> Iterator[Path]:
#         for fqfile in self.files:
#             yield fqfile

#     def __str__(self):
#         return f"FastqArtifact({self.r1}\n{self.r2}\n)"


# @dataclass(frozen=True)
# class TableArtifact(Artifact):
#     path: Path
#     index_column: str | None = None
#     # column_schema: ColumnSchema = field(default_factory=lambda: ColumnSchema({})) # noqa: E501

#     ############################################################################
#     # Artifact protocol, required by Element
#     ############################################################################
#     def resolve(self) -> Path:
#         return self.path

#     def signature(self) -> str:
#         return file_signature(self.path)

#     ############################################################################
#     # Loading
#     ############################################################################
#     @cached_property
#     def stem(self) -> str:
#         return self.path.stem

#     @cached_property
#     def frame(self) -> DataFrame:
#         kwargs = (
#             {"index_col": self.index_column}
#             if self.path.suffix in [".tsv", ".csv"]
#             else {}
#         )
#         frame = read_frame(self.path, **kwargs)
#         if self.index_column and self.index_column in frame.columns:
#             frame.set_index(self.index_column, inplace=True)
#             # frame.index = Index(frame[self.index_column].copy())
#         return frame

#     @cached_property
#     def schema(self) -> ColumnSchema:
#         # return self.column_schema
#         return ColumnSchema.derive_from_columns(
#             self.raw_columns(), index_column=self.index_column
#         )

#     def raw_columns(self) -> list[str]:
#         """Column names on disk, without loading the full table."""
#         raw_columns = [
#             col for col in read_schema(self.path) if col != self.index_column
#         ]
#         return raw_columns

#     ############################################################################
#     # Selection via ColumnSchema
#     ############################################################################

#     def view(
#         self,
#         view: View | None = None,
#     ) -> DataFrame:
#         cols = self.schema.select(view)
#         df = self.frame
#         selected_columns = [col for col in cols if col in df.columns]
#         # difference = set(cols) - set(df.columns)
#         # if difference:
#         #     raise KeyError(
#         #         f"Columns {difference} are in the schema but not in the DataFrame for {self.path}."  # noqa: E501
#         #     )
#         if not cols and view:
#             print(self.frame.head())
#             raise KeyError(
#                 f"No columns match metric={view.metric!r}, sample_id={view.sample_id!r}, pipeline={view.pipeline!r}."  # noqa: E501
#                 f"in {self.path}."
#             )
#         return self.frame[selected_columns]  # type: ignore[return-value]

#     def deselect(
#         self,
#         *,
#         view: View | None = None,
#         keep_untagged: bool = True,
#     ) -> DataFrame:
#         cols = self.schema.deselect(view)
#         if keep_untagged:
#             untagged = [
#                 c for c in self.frame.columns if c not in self.schema
#             ]  # noqa: E501
#             cols = untagged + [c for c in cols if c not in untagged]
#         if not cols and view is not None:
#             raise KeyError(
#                 f"No columns remain after excluding metrics={view.metric!r}, sample_id={view.sample_id!r}, pipeline={view.pipeline!r} in {self.path}."  # noqa: E501
#             )
#         # preserve original column order
#         ordered = [c for c in self.frame.columns if c in cols]
#         return self.frame[ordered]

#     def roles_present(self) -> set[tuple[str, ...] | str]:
#         return {tag.metric for tag in self.schema.values()}

#     ############################################################################
#     # Convenience
#     ############################################################################

#     def with_schema(self, schema: ColumnSchema) -> "TableArtifact":
#         """Return a new TableArtifact pointing at the same path/index
#         but with a different (e.g. merged) ColumnSchema."""
#         return TableArtifact(
#             path=self.path,
#             index_column=self.index_column,
#             # schema=schema,
#         )


# @dataclass(frozen=True)
# class FileArtifact(Artifact):
#     path: Path

#     def resolve(self) -> Path:
#         return self.path

#     def signature(self) -> str:
#         return file_signature(self.path)


# @dataclass(frozen=True)
# class ValueArtifact(Artifact):
#     value: Any

#     def resolve(self) -> Any:
#         return self.value

#     def signature(self) -> str:
#         return stable_hash(self.value)


# @dataclass(frozen=True)
# class DynamicArtifact(Artifact):
#     value: DynamicValue

#     def resolve(self) -> Any:
#         return self.value.resolve()

#     def signature(self) -> str:
#         return self.value.signature


# @dataclass(frozen=True)
# class TransientArtifact(Artifact):
#     value: Any

#     def resolve(self) -> Any:
#         return self.value

#     def signature(self) -> str:
#         return "transient"  # transients are not considered for signature


@dataclass(frozen=True)
class FileSpec:
    "Short specification for an additional output file"

    stem: str
    ext: str


@dataclass(frozen=True)
class PartialOutSpec(Overlay["PartialOutSpec"]):
    """
    Partial specification for the output of an Element, allowing for optional
    overrides.
    """

    stem: str | None = None
    outdir: Path | None = None
    prefix: str | None = None
    suffix: str | None = None
    ext: str | None = None
    exts: tuple[str, ...] | None = None
    additional_output: dict[str, FileSpec] | None = None


@dataclass(frozen=True)
class OutSpec(Overlay[PartialOutSpec]):
    """
    A specification for the output of an Element, including the output
    directory, file name, and extension. Used to generate default names and
    paths for output files and to override defaults specifically.
    """

    stem: str | None = None
    outdir: Path | None = None
    prefixes: list[str] = field(default_factory=list)
    suffixes: list[str] = field(default_factory=list)
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

    def __post_init__(self):
        self.stem = self.stem or "output"
        self.outdir = self.outdir or Path("results")
        output_groups = dict[str, list[Path]] = {}
        for ext in [self.ext] + list(self.exts or ()):
            group_paths = []
            for prefix in self.prefixes:
                for suffix in self.suffixes:
                    fullname = self._make_file_name(
                        prefix=prefix, stem=self.stem, suffix=suffix, ext=ext
                    )
                    fullpath = self._make_path(fullname)
                    group_paths.append(fullpath)
            output_groups[ext] = group_paths
        for key, file_spec in self.additional_output.items():
            fullname = self._make_file_name(
                prefix="", stem=file_spec.stem, suffix="", ext=file_spec.ext
            )
            fullpath = self._make_path(fullname)
            output_groups[key] = [fullpath]
        self.output_groups = output_groups
        missing = [f for f in self._REQUIRED if getattr(self, f) is None]
        if missing:
            raise ValueError(
                f"OutSpec: required field(s) must not be None: {', '.join(missing)}"  # noqa: E501
            )

    def _make_path(self, filename: str) -> Path:
        if self.outdir is None:
            return Path(filename)

        return self.outdir / filename

    def _make_file_name(self, prefix: str, stem: str, suffix: str, ext: str) -> str:
        return f"{prefix}{stem}{suffix}.{ext}"

    def add_output(self, key: str, file_spec: FileSpec) -> OutSpec:
        """
        Add an additional output file specification to the OutSpec.

        Parameters
        ----------
        key : str
            The key for the additional output file.
        file_spec : FileSpec
            The specification for the additional output file.

        Returns
        -------
        OutSpec
            A new OutSpec instance with the additional output file added.
        """
        new_additional_output = dict(self.additional_output)
        new_additional_output[key] = file_spec
        return replace(self, additional_output=new_additional_output)

    def resolve(self, *patches: OutSpec | PartialOutSpec | None) -> OutSpec:
        """
        Resolve the OutSpec by applying the given patches.

        Each patch can be an OutSpec, PartialOutSpec, or None. The patches are
        applied in order, with later patches overriding earlier ones.

        Returns
        -------
        OutSpec
            The resolved OutSpec.
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

    def iter_grouped(self) -> Iterator[tuple[str, list[Path]]]:
        yield (self.ext, self.output_groups[self])

        for key, paths in self.output_groups.items():
            if key == self.ext:
                continue
            for path in paths:
                yield key, path

    def iterfiles(self) -> Iterator[Path]:
        if self.stem is None:
            return iter(())

        # primary
        for fullpath in self.output_groups.get(self.ext, []):
            yield fullpath

        for key in self.output_groups:
            if key == self.ext:
                continue
            for fullpath in self.output_groups[key]:
                yield fullpath

    def filename(self) -> str:
        paths = self.output_groups[self.ext]
        if len(paths) == 0:
            raise ValueError("No output files defined for the primary extension.")
        elif len(paths) > 1:
            raise ValueError("Multiple output files defined for the primary extension.")
        else:
            return paths[0].name

    def path(self) -> Path:
        filename = self.filename()
        return self._make_path(filename)

    @classmethod
    def from_tag(cls, tag: ElementTag, patch: PartialOutSpec | None = None) -> OutSpec:
        """
        Create an OutSpec from an ElementTag, optionally applying a patch.

        Parameters
        ----------
        tag : ElementTag
            The ElementTag to create the OutSpec from.
        patch : PartialOutSpec | None, optional
            An optional patch to apply to the OutSpec, by default None.

        Returns
        -------
        OutSpec
            The created OutSpec.
        """
        base = OutSpec(
            stem=tag.default_name,
            outdir=Path("results") / tag.root,
            ext="parquet",
        )
        return base.patch(patch)
