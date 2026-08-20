from __future__ import annotations

from abc import ABC
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from types import MappingProxyType
from typing import Any, Self

from pandas import DataFrame  # type: ignore[import]

from mmalignments.services.dependencies import (
    Signifiable,
    combined_signature,
    file_signature,
    stable_hash,
)
from mmalignments.services.io import read_frame  # , read_schema


class Artifact(ABC):
    """
    Base class for artifacts.

    Artifacts are immutable objects that represent data or files produced by an
    Element. They can be used to compute signatures and track dependencies.
    """

    pass


class FileArtifact(Artifact, Signifiable):

    @cached_property
    def files(self) -> tuple[Path, ...]:
        "Resolved file paths as a tuple."
        return tuple(ffile for ffile in self)

    @cached_property
    def signature(self) -> str:
        "The combined signature of the files."
        return stable_hash(tuple(file_signature(f) for f in self.files))

    def __iter__(self) -> Iterator[Path]:
        "iterate over resolved file paths"
        for fqfile in self.files:
            yield fqfile


@dataclass(frozen=True)
class FastqArtifact(FileArtifact):
    """
    A FastqArtifact represents a single-end or paired-end FASTQ file(s).
    """

    r1: Path
    r2: Path | None = None
    FASTQ_SUFFIXES: tuple[str, ...] = ("fq", "fastq", "fq.gz", "fastq.gz")

    @cached_property
    def files(self) -> tuple[Path, ...]:
        if self.paired:
            return (self.r1.resolve(), self.r2.resolve())  # type: ignore[return-value]
        else:
            return (self.r1.resolve(),)

    @cached_property
    def stem(self) -> str:
        stem = self.r1.name
        for suffix in self.FASTQ_SUFFIXES:
            rem = "." + suffix
            if stem.endswith(rem):
                stem = stem.removesuffix(rem)
                break
        return stem

    @property
    def paired(self) -> bool:
        return self.r2 is not None

    def __str__(self):
        "String representation of the FastqArtifact."
        return f"FastqArtifact({self.r1}\n{self.r2}\n)"


@dataclass(frozen=True)
class TableArtifact(FileArtifact):
    """
    A TableArtifact represents a tabular data file, such as a CSV or TSV file.
    """

    path: Path
    index_column: str | None = None
    # column_schema: ColumnSchema = field(default_factory=lambda: ColumnSchema({})) # noqa: E501

    ############################################################################
    # Artifact protocol, required by Element
    ############################################################################

    @cached_property
    def signature(self) -> str:
        return file_signature(self.path)

    @cached_property
    def files(self) -> tuple[Path, ...]:
        return (self.path.resolve(),)

    @cached_property
    def stem(self) -> str:
        return self.path.stem

    ############################################################################
    # Loading
    ############################################################################
    def resolve(self) -> Path:
        return self.path.resolve()

    @cached_property
    def frame(self) -> DataFrame:
        return read_frame(self.resolve())

    def read(
        self, drop_unnamed_columns: bool = True, index_column: str | None = None
    ) -> DataFrame:
        kwargs = {}
        if index_column and self.resolve().suffix in [".tsv", ".csv"]:
            kwargs = {"index_col": index_column}
        frame = read_frame(
            self.resolve(), drop_unnamed_columns=drop_unnamed_columns, **kwargs
        )
        if index_column and index_column in frame.columns:
            frame.set_index(index_column, inplace=True)
            # frame.index = Index(frame[index_column].copy())
        return frame


class ArtifactSet(Mapping[str, Any], Signifiable):
    """
    Immutable collection of artifacts with one distinguished primary artifact.

    The primary artifact is always accessible via the key ``"primary"``.
    Additionally, it can be accessed via ``primary_name``.

    Examples
    --------
    >>> artifacts = ArtifactSet(Path("sample.bam"))
    >>> artifacts["primary"]
    PosixPath("sample.bam")

    >>> artifacts["bam"]
    PosixPath("sample.bam")
    """

    def __init__(
        self,
        primary: Any,
        /,
        primary_name: str | None = None,
        **extras: Any,
    ) -> None:
        """
        Create an ArtifactSet.

        Parameters
        ----------
        primary:
            The primary artifact.

        primary_name:
            Optional name for the primary artifact.
            If omitted, the name is inferred from the artifact.

        extras:
            Additional named artifacts.
        """

        if "primary" in extras:
            raise KeyError("'primary' is reserved for the primary artifact.")

        resolved_name = primary_name or self._infer_primary_name(primary)

        if resolved_name in extras:
            raise KeyError(
                f"Primary name '{resolved_name}' collides with extra artifact."
            )

        self._primary = primary
        self._primary_name = resolved_name
        self._extras = MappingProxyType(dict(extras))

    ############################################################################
    # Properties
    ############################################################################

    @property
    def primary(self) -> Any:
        """
        Return the primary artifact.
        """

        return self._primary

    @property
    def primary_name(self) -> str:
        """
        Return the canonical name of the primary artifact.
        """

        return self._primary_name

    @property
    def extras(self) -> Mapping[str, Any]:
        """
        Return the additional artifacts.
        """

        return self._extras

    @cached_property
    def signature(self) -> str:
        "The signature used to invalidate this Element and triggerr a re-run."
        return self.calculate_signature()

    ############################################################################
    # Signature computation
    ############################################################################

    @cached_property
    def signature_data(self) -> Mapping[str, Any]:
        "The data used to compute the signature of this Element."
        return dict(sorted(self.items()))

    def calculate_signature(self) -> str:
        return combined_signature(*self.signature_data.values())

    ############################################################################
    # Convenience methods for file-like artifacts
    ############################################################################

    @cached_property
    def files(self) -> tuple[Path, ...]:
        """
        Return resolved output paths of file-like artifacts.
        """
        all_files = ()
        for artifact in self.values():
            if hasattr(artifact, "files"):
                all_files += artifact.files
            elif isinstance(artifact, Path):
                all_files += (artifact.resolve(),)
        return all_files

    def iterfiles(self) -> Iterator[Path]:
        """
        Iterate over all artifact files.
        """
        for filepath in self.files:
            yield filepath

    ############################################################################
    # ArtifactSet manipulation and creation
    ############################################################################

    @staticmethod
    def _infer_primary_name(value: Any) -> str:
        """
        Infer a useful name for the primary artifact.

        For file-like artifacts the file suffix is used.

        Examples
        --------
        ``sample.bam`` -> ``"bam"``
        ``sample.parquet`` -> ``"parquet"``

        If no meaningful name can be inferred, ``"primary"`` is returned.
        """

        path: Path | None = None

        if isinstance(value, Path):
            path = value

        elif hasattr(value, "resolve"):
            resolved = value.resolve()
            if isinstance(resolved, Path):
                path = resolved

        if path is not None:
            suffix = path.suffix.removeprefix(".")

            if suffix:
                return suffix

        return "primary"

    @classmethod
    def from_any(
        cls,
        value: ArtifactSet | Mapping[str, Any],
        primary_key: str = "primary",
    ) -> ArtifactSet:
        """
        Convert a mapping into an ArtifactSet.

        If no explicit primary key exists, common file formats are checked
        before falling back to the first item.
        """

        if isinstance(value, ArtifactSet):
            return value

        if not isinstance(value, Mapping):
            raise TypeError(f"Cannot create ArtifactSet from {type(value)!r}")

        if primary_key in value:
            primary = value[primary_key]
            primary_name = primary_key

        else:
            for candidate in ("bam", "parquet", "tsv", "fastq"):
                if candidate in value:
                    primary_name = candidate
                    primary = value[candidate]
                    break
            else:
                primary_name, primary = next(iter(value.items()))

        return cls(
            primary,
            primary_name=primary_name,
            **{
                key: artifact
                for key, artifact in value.items()
                if key != primary_name  # noqa: E501
            },
        )

    def with_extra(self, key: str, value: Any) -> Self:
        """
        Return a new ArtifactSet with one additional artifact.

        Raises
        ------
        KeyError
            If the artifact name already exists.
        """

        if key in self:
            raise KeyError(f"Artifact '{key}' already exists.")

        return type(self)(
            self._primary,
            primary_name=self._primary_name,
            **self._extras,
            **{key: value},
        )

    def with_extras(self, extras: Mapping[str, Any]) -> Self:
        """
        Return a new ArtifactSet with additional artifacts.

        Raises
        ------
        KeyError
            If any artifact name already exists.
        """

        conflicts = set(extras) & set(self)

        if conflicts:
            raise KeyError(f"Artifacts already exist: {sorted(conflicts)}")

        return type(self)(
            self._primary,
            primary_name=self._primary_name,
            **self._extras,
            **dict(extras),
        )

    def merge(self, other: ArtifactSet | Mapping[str, Any]) -> Self:
        """
        Merge another ArtifactSet into this one.

        The primary artifact of ``self`` is preserved.
        The primary artifact of ``other`` becomes an extra artifact if it has
        another name.

        Raises
        ------
        KeyError
            If artifact names collide.
        """

        if not isinstance(other, ArtifactSet):
            other = ArtifactSet.from_any(other)

        extras = dict(self._extras)

        other_primary_name = other.primary_name

        if other_primary_name != self.primary_name:
            if other_primary_name in extras:
                raise KeyError(
                    f"Artifact '{other_primary_name}' already exists."
                )  # noqa: E501

            extras[other_primary_name] = other.primary

        conflicts = set(extras) & set(other.extras)

        if conflicts:
            raise KeyError(f"Artifacts already exist: {sorted(conflicts)}")

        extras.update(other.extras)

        return type(self)(
            self.primary,
            primary_name=self.primary_name,
            **extras,
        )

    ############################################################################
    # Mapping interface
    ############################################################################

    def __getitem__(self, key: str) -> Any:
        """
        Return an artifact by name.

        ``"primary"`` always resolves to the primary artifact.
        The primary can additionally be accessed by ``primary_name``.
        """

        if key == "primary" or key == self._primary_name:
            return self._primary

        return self._extras[key]

    def __iter__(self) -> Iterator[str]:
        """
        Iterate over artifact names.
        """

        yield self._primary_name
        yield from self._extras

    def __len__(self) -> int:
        """
        Return the number of named artifacts.

        The primary artifact counts once, even though it has two possible
        aliases.
        """

        return 1 + len(self._extras)

    def __contains__(self, key: object) -> bool:
        """
        Check whether an artifact name exists.
        """

        return (
            key == "primary" or key == self._primary_name or key in self._extras
        )  # line-too-long

    def __repr__(self) -> str:
        """
        Return a readable representation.
        """

        extras = ", ".join(
            f"{key}={value!r}" for key, value in self._extras.items()
        )  # extra string representation

        primary = f"{self._primary_name}={self._primary!r}"

        if extras:
            return f"ArtifactSet({primary}, {extras})"

        return f"ArtifactSet({primary})"
