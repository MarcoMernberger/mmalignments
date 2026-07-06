from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Callable, Iterable, Literal

_GENE_ANNOTATIONS = {
    "gene_stable_id": "META",
    "name": "META",
    "chr": "META",
    "start": "META",
    "stop": "META",
    "strand": "META",
    "tss": "META",
    "tes": "META",
    "biotype": "META",
    "transcript_stable_ids": "META",
}


@dataclass(frozen=True)
class View:
    """
    metric : str | tuple[str, ...] | None, optional
        The metric prefix to filter by. If None, do not filter by metric.
    sample_id : str | tuple[str, ...] | None, optional
        The sample ID or sequence of sample IDs to filter by. If None, do
        not filter by sample ID.
    """

    # role: Role | None = None
    metric: str | tuple[str, ...] | list[str] | None = None
    sample_id: str | tuple[str, ...] | None = None
    pipeline: str | None = None  #  e.g. "STAR_UMI-tools_dedup-directional"


_VIEWS = {  # often used, for convenience
    "gu": View(metric=("Count", "GU")),
    "gud": View(metric=("Count", "GU", "Dedup")),
    "esc": View(metric=("Count", "ES")),
    "escd": View(metric=("Count", "ES", "Dedup")),
    "vst": View(metric=("VST",)),
}


KNOWN_METRICS: dict[str, str] = {
    "gu": "Gene unstranded count",
    "gud": "Gene unstranded count deduplicated",
    "esc": "Exon stranded count",
    "escd": "Exon stranded count, deduplicated",
    "vst": "Variance-stabilised transformation (DESeq2)",
    "zscore": "Z-standardized score",
}


def describe_metric(segment: str) -> str:
    return KNOWN_METRICS.get(segment, f"(unknown metric segment: {segment!r})")


_CATEGORY_METRIC = {
    ("exon, protein coding, stranded", True): _VIEWS["escd"].metric,  # dedup
    ("exon, protein coding, stranded", False): _VIEWS["esc"].metric,  # no dedup
    ("gene unstranded", True): _VIEWS["gud"].metric,  # dedup
    ("gene unstranded", False): _VIEWS["gu"].metric,  # no dedup
}

################################################################################
# Legacy column name parsing
################################################################################


def _legacy_to_tag(colname: str, sample_names: list[str]) -> ColumnTag | None:
    if colname in _GENE_ANNOTATIONS:
        return ColumnTag(
            metric=(colname,),
            sample_id=None,  # role=_GENE_ROLES[colname]
            pipeline=_GENE_ANNOTATIONS[colname],  # e.g. "META"
        )  # leave as is
    else:
        m = re.match(r"^(?P<category>.+?)\s+tag count\s+(?P<rest>.+)$", colname)
        if not m:
            return None
        category = m.group("category").strip().lower()
        rest = m.group("rest").strip()

        # längste bekannte sample_id matchen, die 'rest' als Prefix hat
        sample_id = max(
            (s for s in sample_names if rest == s or rest.startswith(s + "_")),
            key=len,
            default=None,
        )
        if sample_id is None:
            raise ValueError(
                f"Could not match any known sample_id in column: {colname!r}"
            )

        pipeline = rest[len(sample_id) :].lstrip(
            "_"
        )  # e.g. "STAR_UMI-tools_dedup-directional"

        is_dedup = "dedup" in pipeline.lower()
        role_metric = _CATEGORY_METRIC.get((category, is_dedup))
        if role_metric is None:
            raise ValueError(
                f"Unknown category/dedup combination: {(category, is_dedup)}"
            )
        return ColumnTag(
            metric=role_metric,
            sample_id=sample_id,
            # role=role_metric[0],
            pipeline=None,  # pipeline
        )


################################################################################
# use Column Tags for selecting
################################################################################

# A regular expression pattern to match canonical column names in the format:
# "<metric> (<sample>)[pipeline]"
_COL_PATTERN = re.compile(
    r"""
    ^
    (?P<metric>[^\[(]+?)          # metric
    (?:\s*\((?P<sample>[^()]+)\))? # optional (sample)
    (?:\s*\[(?P<pipeline>[^\[\]]+)\])? # optional [pipeline]
    $
    """,
    re.VERBOSE,
)


@dataclass(frozen=True)
class ColumnTag:
    """
    A tag that describes a column in a DataFrame, including its metric,
    sample ID, and pipeline. Used for easy selection and filtering of columns
    based on their characteristics.
    """

    metric: tuple[str, ...]
    sample_id: str | None = None
    pipeline: str | None = None  # "STAR_UMI-tools_dedup-directional", optional

    def encode(self) -> str:
        """
        Generate a canonical column name from the tag's attributes.

        Returns
        -------
        str
            The canonical column name generated from the tag's attributes.
        """
        column_name = f"{'.'.join(self.metric)}"
        if self.sample_id:
            column_name = f"{column_name} ({self.sample_id})"
        if self.pipeline:
            column_name = f"{column_name} [{self.pipeline}]"
        return column_name

    @classmethod
    def decode(cls, colname: str) -> ColumnTag | None:
        """
        Decode a canonical column name into a ColumnTag.

        The column name must follow the format generated by ColumnTag.encode().

        Parameters
        ----------
        colname : str
            The canonical column name to decode.
        # role : Role
        #     The role associated with the column.

        Returns
        -------
        ColumnTag | None
            The decoded ColumnTag if the column name matches the expected
            format, otherwise None.
        """
        m = _COL_PATTERN.match(colname.strip())
        if not m:
            return None

        sample = m.group("sample")
        pipeline = m.group("pipeline")

        return cls(
            metric=tuple(
                part.strip() for part in m.group("metric").strip().split(".")
            ),  # noqa: E501
            sample_id=sample.strip() if sample else None,
            pipeline=pipeline.strip() if pipeline else None,
        )


class ColumnSchema(Mapping[str, ColumnTag]):
    """
    A schema that describes the columns in a DataFrame using ColumnTags.
    Provides methods for selecting and deselecting columns based on their tags.
    """

    __slots__ = ("_tags",)

    def __init__(self, tags: Mapping[str, ColumnTag], index: str | None = None):
        """
        Initialize the ColumnSchema with a mapping of columnnames to ColumnTags.

        Parameters
        ----------
        tags : Mapping[str, ColumnTag]
            A mapping from column names to their corresponding ColumnTags.
        """
        object.__setattr__(self, "_tags", MappingProxyType(dict(tags)))
        self.index_column = index

    def __getitem__(self, key: str) -> ColumnTag:
        return self._tags[key]

    def __iter__(self):
        return iter(self._tags)

    def __len__(self) -> int:
        return len(self._tags)

    ############################################################################
    # Selection and filtering methods
    ############################################################################

    def with_columns(self, tags: Mapping[str, ColumnTag]) -> ColumnSchema:
        overlap = set(tags) & set(self._tags)
        if overlap:
            raise KeyError(f"Columns already in schema: {overlap}")
        return ColumnSchema({**self._tags, **tags})

    def _matches(
        self,
        tag: ColumnTag,
        *,
        # role: Role | None,
        metric_prefix: tuple[str, ...] | None,
        samples: set[str] | None,
    ) -> bool:
        """
        Check if a column tag matches the given filtering criteria.

        This method is used internally by the select and deselect methods to
        determine if a column should be included based on its role, metric
        prefix and sample ID.

        Parameters
        ----------
        tag : ColumnTag
            The ColumnTag to check against the filtering criteria.
        role : Role | None
            The role to filter by. If None, do not filter by role.
        metric_prefix : tuple[str, ...] | None
            The metric prefix to filter by. If None, do not filter by metric.
        samples : set[str] | None
            The set of sample IDs to filter by. If None, do not filter by
            sample ID.

        Returns
        -------
        bool
            True if the tag matches all specified criteria, False otherwise.
        """
        # if role is not None and tag.role != role:
        #     return False
        if (
            metric_prefix is not None
            and tag.metric[: len(metric_prefix)] != metric_prefix
        ):
            return False
        if samples is not None and tag.sample_id not in samples:
            return False
        return True

    def select(
        self,
        view: View | None = None,
    ) -> list[str]:
        """
        Select columns based on their tags.

        This method allows filtering columns by their role, metric prefix, and
        sample ID.

        Parameters
        ----------

        Returns
        -------
        list[str]
            A list of column names that match the specified filtering criteria.
        """
        if view is None:
            print("view is None, returning all columns")
            return list(self._tags.keys())
        metric_prefix = (
            (view.metric,) if isinstance(view.metric, str) else view.metric
        )  # noqa: E501
        samples = (
            {view.sample_id}
            if isinstance(view.sample_id, str)
            else (set(view.sample_id) if view.sample_id is not None else None)
        )
        for col, tag in self._tags.items():
            print(
                f"Checking column: {col}, tag: {tag}: matches",
                self._matches(tag, metric_prefix=metric_prefix, samples=samples),
            )
        return [
            col
            for col, tag in self._tags.items()
            if self._matches(
                tag,
                metric_prefix=metric_prefix,
                samples=samples,
            )
        ]

    def deselect(
        self,
        view: View | None = None,
    ) -> list[str]:
        """
        Deselect columns based on their tags.

        This method allows filtering columns by their role, metric prefix, and
        sample ID, and returns the columns that do NOT match the specified
        criteria.


        Parameters
        ----------
        view : View | None, optional
            The view containing filtering criteria. If None, do not filter by
            role, metric, or sample ID.

        Returns
        -------
        list[str]
            A list of column names that do NOT match the specified filtering
            criteria.
        """
        if view is None:
            return []
        metric_prefix = (
            (view.metric,) if isinstance(view.metric, str) else view.metric
        )  # noqa: E501
        samples = (
            {view.sample_id}
            if isinstance(view.sample_id, str)
            else (set(view.sample_id) if view.sample_id is not None else None)
        )
        return [
            col
            for col, tag in self._tags.items()
            if not self._matches(
                tag,
                metric_prefix=metric_prefix,
                samples=samples,
            )
        ]

    def generate_new_tag(
        self, column_name, old_tag, add, append, new_pipeline
    ) -> ColumnTag:
        new_tag = ColumnTag(
            metric=old_tag.metric + add if append else add,
            sample_id=old_tag.sample_id,
            pipeline=new_pipeline,
        )
        return new_tag

    def derive(
        self,
        columns: Iterable[str],
        *,
        view_new: View | None = None,
        append: bool = False,
    ) -> tuple[ColumnSchema, dict[str, str]]:
        """
        Derive new columns from existing ones by modifying their tags.
        Also returns a mapping from old column names to new column names, to
        rename newly created columns.

        Parameters
        ----------
        columns : Iterable[str]
            Existing column names that must be tagged in the schema.
        add_metric : str | tuple[str, ...]
            Appended to tag.metric, e.g., "vst" -> metric + ("vst",)
        # role : Role
        #     New role for the derived columns.
        pipeline : str | None
            If not specified (sentinel), the pipeline field from the original
            tag is used. Explicit None deletes it, a string sets it anew.
        append : bool, optional
            If True, the new metric is appended to the existing metric. If
            False the new metric replaces the existing metric. Default is False.

        Returns
        -------
        tuple[ColumnSchema, dict[str, str]]
            The derived column schema and a mapping from old column names to new
            column names.
        """
        add = (
            (view_new.metric,)
            if isinstance(view_new.metric, str)
            else tuple(view_new.metric)
        )

        new_tags: dict[str, ColumnTag] = {}
        rename: dict[str, str] = {}
        for col in columns:
            old_tag = self.get(col)
            if old_tag is None:
                raise KeyError(
                    f"Column {col!r} not found in schema; cannot derive from it."  # noqa: E501
                )

            new_pipeline = (
                view_new.pipeline
                if view_new and view_new.pipeline is not None
                else old_tag.pipeline
            )
            sample_id = (
                view_new.sample_id
                if view_new and view_new.sample_id is not None
                else old_tag.sample_id
            )
            new_tag = ColumnTag(
                metric=old_tag.metric + add if append else add,
                sample_id=sample_id,
                pipeline=new_pipeline,
            )
            new_name = new_tag.encode()
            if new_name in new_tags:
                raise ValueError(
                    f"Duplicate derived column name: {new_name!r} (from {col!r})"  # noqa: E501
                )

            new_tags[new_name] = new_tag
            rename[col] = new_name

        return ColumnSchema(new_tags), rename

    @staticmethod
    def derive_from_columns(
        columns: Iterable[str],
        *,
        view_new: View | None = None,
        append: bool = False,
    ) -> tuple[ColumnSchema, dict[str, str]]:
        """
        Derive new columns from existing ones by modifying their tags.
        Also returns a mapping from old column names to new column names, to
        rename newly created columns.

        Parameters
        ----------
        columns : Iterable[str]
            Existing column names that must be tagged in the schema.
        add_metric : str | tuple[str, ...]
            Appended to tag.metric, e.g., "vst" -> metric + ("vst",)
        # role : Role
        #     New role for the derived columns.
        pipeline : str | None
            If not specified (sentinel), the pipeline field from the original
            tag is used. Explicit None deletes it, a string sets it anew.
        append : bool, optional
            If True, the new metric is appended to the existing metric. If
            False, the new metric replaces the existing metric.
            Default is False.

        Returns
        -------
        tuple[ColumnSchema, dict[str, str]]
            The derived column schema and a mapping from old column names to new
            column names.
        """

        new_tags: dict[str, ColumnTag] = {}

        for col in columns:
            print(col)
            new_tag = ColumnTag.decode(col)
            new_name = new_tag.encode()
            if new_name in new_tags:
                raise ValueError(
                    f"Duplicate derived column name: {new_name!r} (from {col!r})"  # noqa: E501
                )
            new_tags[new_name] = new_tag

        return ColumnSchema(new_tags)

    @staticmethod
    def derive_from_view(view_out: View) -> ColumnSchema:
        """
        Derive a new ColumnSchema based on a future output view.

        Parameters
        ----------
        view_out : View
            The output view containing the new metric and pipeline information.
        append : bool, optional
            If True, the new metric is appended to the existing metric. If
            False, it replaces the existing metric. By default False.

        Returns
        -------
        ColumnSchema
            A new ColumnSchema instance with the updated metric and pipeline.
        """
        new_tags: dict[str, ColumnTag] = {}
        if view_out.sample_id is None:
            return ColumnSchema({})
        metrics_to_produce = [""]
        if isinstance(view_out.metric, list):
            metrics_to_produce = view_out.metric
        elif isinstance(view_out.metric, (str, tuple)):
            metrics_to_produce = [view_out.metric]

        for sample in view_out.sample_id:
            for metric in metrics_to_produce:
                new_tag = ColumnTag(
                    metric=metric if isinstance(metric, tuple) else (metric,),
                    sample_id=sample,
                    pipeline=view_out.pipeline,
                )
                new_name = new_tag.encode()
                if new_name in new_tags:
                    raise ValueError(
                        f"Duplicate derived column name: {new_name!r} (from sample {sample!r})"  # noqa: E501
                    )
                new_tags[new_name] = new_tag
        return ColumnSchema(new_tags)

    def apply_out_view(
        tag: ColumnTag, view_out: View, append: bool = False
    ) -> ColumnTag:
        """
        Create a new ColumnTag based on a future output view.

        Parameters
        ----------
        tag : ColumnTag
            The old Tag to be modified.
        view_out : View
            The output view containing the new metric and pipeline information.
        append : bool, optional
            If True, the new metric is appended to the existing metric. If
            False, it replaces the existing metric. By default False.

        Returns
        -------
        ColumnTag
            A new ColumnTag instance with the updated metric and pipeline.
        """
        add = (
            (view_out.metric,)
            if isinstance(view_out.metric, str)
            else (view_out.metric or ())
        )
        return ColumnTag(
            metric=tag.metric + add,
            sample_id=tag.sample_id,
            pipeline=(
                view_out.pipeline
                if view_out.pipeline is not None
                else tag.pipeline  # noqa: E501
            ),
        )

    def drop(self, columns: Iterable[str]) -> ColumnSchema:
        """
        Drop columns from the schema.

        Returns a new ColumnSchema instance without the specified columns.

        Parameters
        ----------
        columns : Iterable[str]
            Column names to drop from the schema.

        Returns
        -------
        ColumnSchema
            A new ColumnSchema instance without the specified columns.

        Raises
        ------
        KeyError
            If any of the specified columns are not present in the schema.
        """
        to_drop = set(columns)
        missing = to_drop - set(self._tags)
        if missing:
            raise KeyError(f"Cannot drop unknown columns: {missing}")
        return ColumnSchema(
            {k: v for k, v in self._tags.items() if k not in to_drop}
        )  # noqa: E501

    @classmethod
    def build(cls, columns: Iterable[str]) -> "ColumnSchema":
        """
        Build a ColumnSchema for columns that all share the same role.

        Parameters
        ----------
        columns : Iterable[str]
            column names from DataFrame or other source.
        role : Role
            The role to assign to all columns in the schema.

        Returns
        -------
        ColumnSchema
            New ColumnSchema instance with the given columns and role.
        """
        tags = {}
        for col in columns:
            tag = ColumnTag.decode(col)  # role=role
            if tag is not None:
                tags[col] = tag
        return cls(tags)

    def merge(
        self,
        other: "ColumnSchema",
        *,
        on_conflict: Literal["raise", "keep_self", "keep_other"] = "raise",
    ) -> "ColumnSchema":
        """
        Merge this ColumnSchema with another, handling conflicts according to
        the specified strategy.

        Parameters
        ----------
        other : ColumnSchema
            The other ColumnSchema to merge with.
        on_conflict : Literal["raise", "keep_self", "keep_other"], optional
            Strategy to handle conflicts:
            - "raise": raise an error if there are conflicting columns.
            - "keep_self": keep columns from this schema in case of conflict.
            - "keep_other": keep columns from the other schema in case of
               conflict.
            By default "raise".

        Returns
        -------
        ColumnSchema
            A new ColumnSchema instance that is the result of merging this
            schema with the other.

        Raises
        ------
        KeyError
            If there are conflicting columns and on_conflict is set to "raise".
        """
        overlap = set(self._tags) & set(other._tags)
        if overlap and on_conflict == "raise":
            raise KeyError(f"Column tag conflict on merge: {sorted(overlap)}")
        if on_conflict == "keep_self":
            merged = {**other._tags, **self._tags}
        else:  # keep_other oder kein overlap
            merged = {**self._tags, **other._tags}
        return ColumnSchema(merged)

    @classmethod
    def from_legacy(
        cls,
        columns: Iterable[str],
        sample_names: list[str],
        parser: Callable[[str, list[str]], ColumnTag | None] = _legacy_to_tag,
    ) -> tuple["ColumnSchema", dict[str, str]]:
        """
        Convert legacy column names to a new ColumnSchema and provide a rename
        mapping.

        This method parses legacy column names using the provided parser,
        constructs a new ColumnSchema with the recognized columns, and returns
        a mapping from the original column names to the new canonical names.

        Returns
        -------
        tuple[ColumnSchema, dict[str, str]]
            A tuple containing the new ColumnSchema and a rename mapping from
            old column names to new canonical names.

        Raises
        ------
        ValueError
            If there are duplicate columns after renaming.
        """
        tags: dict[str, ColumnTag] = {}
        rename: dict[str, str] = {}
        for col in columns:
            tag = parser(col, sample_names)
            if tag is None:
                continue
            new_name = tag.encode()
            if new_name in tags:
                raise ValueError(
                    f"Duplicate column after rename: {new_name!r} (from {col!r})"  # noqa: E501
                )
            tags[new_name] = tag
            rename[col] = new_name
        return cls(tags), rename
