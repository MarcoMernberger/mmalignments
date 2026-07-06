from __future__ import annotations

from abc import ABC
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from types import MappingProxyType
from typing import Any, Iterator, Literal, Mapping

from pandas import DataFrame

from mmalignments.core.annotations import ColumnSchema, View
from mmalignments.models.tags import ElementTag
from mmalignments.services.dependencies import (
    DynamicValue,
    file_signature,
    stable_hash,
)
from mmalignments.services.io import read_frame, read_schema

###############################################################################
# External Wrapper
###############################################################################

ArtifactType = Literal["file", "value"]
ArtifactLifeTime = Literal["persistent", "transient", "ephemeral"]


@dataclass(frozen=True)
class OutputSpec:
    filename: str | None = None
    outdir: Path | None = None
    additional_extensions: list[str] | None = None
    ext: str | None = None  # e.g. parquet


class ArtifactSet(Mapping):
    __slots__ = ("_primary", "_primary_name", "_extras")

    def __init__(
        self, primary: Any, /, primary_name: str | None = None, **extras: Any
    ):  # noqa: E501
        if "primary" in extras:
            raise KeyError("'primary' is reserved for the primary artifact.")
        if primary_name is not None and primary_name in extras:
            raise KeyError(
                f"'{primary_name}' collides with an extra artifact key."
            )  # noqa: E501
        object.__setattr__(self, "_primary", primary)
        object.__setattr__(self, "_extras", MappingProxyType(dict(extras)))
        object.__setattr__(self, "_primary_name", primary_name)

    def __getitem__(self, key: str) -> Any:
        if key == "primary" or (
            self._primary_name is not None and key == self._primary_name
        ):
            return self._primary
        return self._extras[key]

    def __contains__(self, key: object) -> bool:
        return (
            key == "primary" or key == self._primary_name or key in self._extras
        )  # noqa: E501

    def __iter__(self) -> Iterator[str]:
        yield (
            self._primary_name if self._primary_name is not None else "primary"
        )  # noqa: E501
        yield from self._extras

    def __len__(self) -> int:
        return 1 + len(self._extras)

    @property
    def primary(self) -> Any:
        return self._primary

    @property
    def primary_name(self) -> str | None:
        return self._primary_name

    def with_extra(self, key: str, value: Any) -> ArtifactSet:
        if key == "primary" or key == self._primary_name or key in self._extras:
            raise KeyError(f"Artifact '{key}' already exists.")
        return ArtifactSet(
            self._primary,
            primary_name=self._primary_name,
            **self._extras,
            **{key: value},
        )

    def __repr__(self) -> str:
        extras = ", ".join(f"{k}={v!r}" for k, v in self._extras.items())
        name = self._primary_name or "primary"
        return f"ArtifactSet({name}={self._primary!r}{', ' + extras if extras else ''})"  # noqa: E501

    @classmethod
    def from_any(
        cls, x: "ArtifactSet | Mapping[str, Any]", primary_key: str = "primary"
    ) -> "ArtifactSet":
        """
        Convert a mapping to an ArtifactSet for backwards compatibility. If
        the mapping does not contain a key for the primary artifact, a fallback
        is attempted using the first available value or a known extension.
        """
        if isinstance(x, ArtifactSet):
            return x
        if isinstance(x, Mapping):
            if primary_key not in x:
                if "parquet" in x:
                    primary_name = "parquet"
                    primary = x[primary_name]
                elif "tsv" in x:
                    primary_name = "tsv"
                    primary = x[primary_name]
                else:
                    primary_name, primary = next(iter(x.items()), (None, None))
                if primary is None:
                    raise KeyError(
                        f"Expected key '{primary_key}' in mapping. No fallback found for primary artifact in {list(x.keys())}."  # noqa: E501
                    )
            else:
                primary = x[primary_key]
                primary_name = primary_key

            return cls(
                primary,
                primary_name=primary_name,
                **{k: v for k, v in x.items() if k != primary_name},
            )
        raise TypeError(f"Cannot build ArtifactSet from {type(x)!r}")

    @classmethod
    def generate_file_artifacts(
        cls,
        tag: ElementTag,
        # column_schema: ColumnSchema,
        infile: Path | None,
        spec: OutputSpec | None = None,
        default_dir: Path | None = None,
        index_column: str | None = None,
    ) -> tuple[ArtifactSet, list[Path]]:
        default_dir = default_dir or Path("results")
        spec = spec or OutputSpec()
        out_dir = Path(spec.outdir or infile.parent if infile else default_dir)
        filename = spec.filename or tag.default_name
        if spec.ext:
            filename = f"{filename}.{spec.ext}"
        output_file = out_dir / filename
        extensions = spec.additional_extensions or []
        additional_artifacts = {}
        for extension in extensions:
            additional_artifacts[extension] = output_file.with_suffix(
                f".{extension}"
            )  # noqa: E501
        artifacts = ArtifactSet(
            TableArtifact(
                output_file, index_column=index_column
            ),  # , column_schema=column_schema
            **additional_artifacts,
        )
        return artifacts, [
            a.resolve() for a in artifacts.values() if hasattr(a, "resolve")
        ]


class Artifact(ABC):
    def resolve(self) -> Any:
        raise NotImplementedError

    def signature(self) -> str:
        raise NotImplementedError


# @dataclass(frozen=True)
# class TableArtifact(Artifact):
#     path: Path
#     index_column: str | None = None
#     schema: dict[str, str] = field(default_factory=dict)
#     _df_cache: DataFrame | None = field(default=None, init=False, repr=False)

#     def resolve(self) -> Path:
#         return self.path

#     def signature(self) -> str:
#         return file_signature(self.path)

#     def load(self) -> DataFrame:
#         return read_frame(self.path)

#     # ------------------------------------------------------------------
#     # column_roles helpers
#     # ------------------------------------------------------------------

#     def columns_for_role(self, role: str) -> list[str]:
#         """Return all column names assigned to *role*."""
#         return [col for col, r in self.column_roles.items() if r == role]

#     def roles_present(self) -> set[str]:
#         """Return the set of distinct roles registered for this table."""
#         return set(self.column_roles.values())

#     # ------------------------------------------------------------------
#     # Views
#     # ------------------------------------------------------------------

#     @cached_property
#     def df(self):
#         df = read_frame(self.path)
#         if self.index_column:
#             df.set_index(self.index_column, inplace=True)
#         return df

#     def view(self, role: str | None = None) -> DataFrame:
#         """Return a :class:`~pandas.DataFrame` with only the columns that
#         have the given *role*.

#         The index of the full table is preserved.

#         Parameters
#         ----------
#         role : str | None
#             Semantic role to filter for, e.g. ``"vst"``, ``"annotation"``.
#             If None, return the full DataFrame.

#         Returns
#         -------
#         pd.DataFrame
#             Subset of :attr:`df` with only the matching columns.

#         Raises
#         ------
#         KeyError
#             If no columns are registered for *role*.
#         """
#         if role is None:
#             return self.df
#         cols = self.columns_for_role(role)
#         if not cols:
#             raise KeyError(
#                 f"No columns with role {role!r} in {self.path}.\n"
#                 f"Available roles: {sorted(self.roles_present())}"
#             )
#         return self.df[cols]

#     def roles(
#         self,
#         *roles: str,
#         extra_columns: Iterable[str] | None = None,
#     ) -> DataFrame:
#         """Return a DataFrame with columns matching any of the given *roles*.

#         Parameters
#         ----------
#         *roles : str
#             One or more role names to include.
#         extra_columns : Iterable[str] | None
#             Additional column names to include regardless of their role.

#         Returns
#         -------
#         pd.DataFrame
#             Combined subset; column order follows the original table.
#         """
#         wanted: set[str] = set()
#         for role in roles:
#             wanted.update(self.columns_for_role(role))
#         if extra_columns:
#             wanted.update(extra_columns)
#         # preserve original column order
#         ordered = [c for c in self.df().columns if c in wanted]
#         return self.df()[ordered]

#     def schema(self) -> list[str]:
#         """Return the list of column names in the table."""
#         return read_schema(self.path)


@dataclass(frozen=True)
class TableArtifact(Artifact):
    path: Path
    index_column: str | None = None
    # column_schema: ColumnSchema = field(default_factory=lambda: ColumnSchema({})) # noqa: E501

    ############################################################################
    # Artifact protocol, required by Element
    ############################################################################
    def resolve(self) -> Path:
        return self.path

    def signature(self) -> str:
        return file_signature(self.path)

    ############################################################################
    # Loading
    ############################################################################

    @cached_property
    def frame(self) -> DataFrame:
        kwargs = (
            {"index_col": self.index_column}
            if self.path.suffix in [".tsv", ".csv"]
            else {}
        )
        frame = read_frame(self.path, **kwargs)
        print(f"Loaded DataFrame from {self.path} with shape {frame.shape}")
        print(frame.head())
        for col in frame.columns:
            print("col in frame", col)
        if self.index_column and self.index_column in frame.columns:
            frame.set_index(self.index_column, inplace=True)
            # frame.index = Index(frame[self.index_column].copy())
        return frame

    @cached_property
    def schema(self) -> ColumnSchema:
        # return self.column_schema
        return ColumnSchema.derive_from_columns(self.raw_columns())

    def raw_columns(self) -> list[str]:
        """Column names on disk, without loading the full table."""
        print("raw_columns schema", read_schema(self.path))
        return read_schema(self.path)

    ############################################################################
    # Selection via ColumnSchema
    ############################################################################

    def view(
        self,
        view: View | None = None,
    ) -> DataFrame:
        cols = self.schema.select(view)
        print(
            "missing in frame but selected by schema",
            [c for c in cols if c not in self.frame.columns],
        )
        if not cols and view:
            print(self.frame.head())
            raise KeyError(
                f"No columns match metric={view.metric!r}, sample_id={view.sample_id!r}, pipeline={view.pipeline!r}."  # noqa: E501
                f"in {self.path}."
            )
        print("in cols")
        for col in cols:
            print("col", col)
        return self.frame[cols]

    def deselect(
        self,
        *,
        view: View | None = None,
        keep_untagged: bool = True,
    ) -> DataFrame:
        cols = self.schema.deselect(view)
        if keep_untagged:
            untagged = [
                c for c in self.frame.columns if c not in self.schema
            ]  # noqa: E501
            cols = untagged + [c for c in cols if c not in untagged]
        if not cols and view is not None:
            raise KeyError(
                f"No columns remain after excluding metrics={view.metric!r}, sample_id={view.sample_id!r}, pipeline={view.pipeline!r} in {self.path}."  # noqa: E501
            )
        # preserve original column order
        ordered = [c for c in self.frame.columns if c in cols]
        return self.frame[ordered]

    def roles_present(self) -> set[tuple[str, ...]]:
        return {tag.metric for tag in self.schema.values()}

    ############################################################################
    # Convenience
    ############################################################################

    def with_schema(self, schema: ColumnSchema) -> "TableArtifact":
        """Return a new TableArtifact pointing at the same path/index
        but with a different (e.g. merged) ColumnSchema."""
        return TableArtifact(
            path=self.path,
            index_column=self.index_column,
            # schema=schema,
        )


@dataclass(frozen=True)
class FileArtifact(Artifact):
    path: Path

    def resolve(self) -> Path:
        return self.path

    def signature(self) -> str:
        return file_signature(self.path)


@dataclass(frozen=True)
class ValueArtifact(Artifact):
    value: Any

    def resolve(self) -> Any:
        return self.value

    def signature(self) -> str:
        return stable_hash(self.value)


@dataclass(frozen=True)
class DynamicArtifact(Artifact):
    value: DynamicValue

    def resolve(self) -> Any:
        return self.value.resolve()

    def signature(self) -> str:
        return self.value.signature


@dataclass(frozen=True)
class TransientArtifact(Artifact):
    value: Any

    def resolve(self) -> Any:
        return self.value

    def signature(self) -> str:
        return "transient"  # transients are not considered for signature
