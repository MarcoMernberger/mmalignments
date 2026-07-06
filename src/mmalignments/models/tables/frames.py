from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property, wraps
from pathlib import Path
from typing import Any, Callable, Iterable, Literal, Mapping, Protocol

import pandas as pd
from pandas import DataFrame

from mmalignments.core.annotations import _VIEWS, ColumnSchema, ColumnTag, View
from mmalignments.models.artifacts import ArtifactSet, OutputSpec, TableArtifact
from mmalignments.models.elements import (
    Element,
    FileElement,
    element,
    generate_element_key_name,
)
from mmalignments.models.parameters import Params
from mmalignments.models.tags import (
    ElementTag,
    Method,
    PartialElementTag,
    State,
    from_prior,
)
from mmalignments.services.dependencies import (
    depends,
    # collect_code_dependency,
    # file_sig,
    function_hash,
    stable_hash,
)
from mmalignments.services.io import read_frame, read_schema, write_frames

################################################################################
# Table sources
################################################################################


class TableSource(Protocol):
    @cached_property
    def frame(self) -> DataFrame: ...

    @cached_property
    def schema(self) -> ColumnSchema: ...

    @property
    def index_column(self) -> str | None: ...

    def resolve(self) -> Path: ...

    def view(self, view: View | None = None) -> DataFrame: ...


@dataclass(frozen=True)
class PathTableSource:
    path: Path
    index_column: str | None = None
    schema_parser: Callable[[Iterable[str]], ColumnSchema] = ColumnSchema.build

    @cached_property
    def frame(self) -> DataFrame:
        frame = read_frame(self.path)
        if self.index_column:
            frame = frame.set_index(self.index_column)
        return frame

    @cached_property
    def schema(self) -> ColumnSchema:
        return self.schema_parser(self.frame.columns)

    def resolve(self) -> Path:
        return self.path

    def view(self, view: View | None = None) -> DataFrame:
        cols = self.schema.select(view) if view else list(self.frame.columns)
        return self.frame[cols]


################################################################################
# Morphs
################################################################################


@dataclass(frozen=True)
class Morph:
    fn: Callable[[DataFrame, Any, ...], DataFrame]
    view_in: View | None = None
    parameter: dict[str, str | float | int | bool] = field(default_factory=dict)
    # view_out: View | None = None
    # create_new: bool = False
    # add: bool = False
    # append_metric: bool = True

    @classmethod
    def from_callable(
        cls, fn: Callable, *, view: View | None = None
    ) -> Morph:  # noqa: E501
        # if not hasattr(fn, "morph_role"):
        #     raise TypeError(f"{fn.__name__} is not decorated with @morph(...)")   # noqa: E501
        morph_view = getattr(fn, "morph_view", view)
        return cls(fn=fn, view_in=morph_view)

    def __str__(self):
        return f"Morph(fn={self.fn.__name__}, view_in={self.view_in}, view_out={self.view_out}, create_new={self.create_new}, add={self.add}, append_metric={self.append_metric})"  # noqa: E501

    def signature(self):
        params_hash = stable_hash(self.parameter)
        return function_hash(self.fn) + params_hash


def morph(
    _fn: Callable | None = None,
    *,
    view_in: View | None = None,
    view_out: View | None = None,
    create_new: bool = False,
    add: bool = False,
    append_metric: bool = True,
):
    """Decorator turning a function into a Morph."""

    def decorator(fn: Callable[..., DataFrame]) -> Morph:
        @wraps(fn)
        def wrapped(*args, **kwargs):
            return fn(*args, **kwargs)

        return Morph(
            fn=wrapped,
            view_in=view_in,
            # view_out=view_out,
            # create_new=create_new,
            # add=add,
            # append_metric=append_metric,
        )

    if _fn is None:
        return decorator

    return decorator(_fn)


@dataclass(frozen=True)
class MultiMorph(Morph):
    fn: Callable[[Mapping[str, DataFrame], Any, ...], DataFrame]
    views: Mapping[str, View] | None = None

    @classmethod
    def from_callable(
        cls, fn: Callable, *, views: Mapping[str, View] | None = None
    ) -> MultiMorph:
        # if not hasattr(fn, "morph_role"):
        #     raise TypeError(f"{fn.__name__} is not decorated with @morph(...)")   # noqa: E501
        morph_views = getattr(fn, "morph_views", views)
        return cls(
            fn=fn, views={"default": morph_views} if morph_views else None
        )  # noqa: E501

    def __str__(self):
        return f"MultiMorph(fn={self.fn.__name__}, views={self.views})"  # noqa: E501


################################################################################
# Filter functions
################################################################################


def filter_to_threshold(
    threshold: int | float = 10,
    how: Literal["any", "all", "atleast"] = "any",
    atleast: int = 3,
    columns: list[str] | None = None,
    view: View | None = _VIEWS["escd"],
    greater_than: bool = True,
) -> Callable[[DataFrame], DataFrame]:

    @morph(params=[threshold, how])
    def filter_func(df: DataFrame) -> DataFrame:
        count_columns = columns if columns is not None else df.columns.to_list()
        passed = (
            df[count_columns].gt(threshold)
            if greater_than
            else df[count_columns].lt(threshold)
        )

        if how == "any":
            mask = passed.any(axis=1)
        elif how == "all":
            mask = passed.all(axis=1)
        elif how == "atleast":
            mask = passed.sum(axis=1) >= atleast
        else:
            raise ValueError(
                f"Invalid 'how' parameter: {how}. Use 'any', 'all', or 'atleast'."  # noqa: E501
            )

        return df[mask]

    return filter_func


def filter_abs_threshold(
    columns: list[str] | str,
    threshold: int | float = 1,
    how: Literal["any", "all", "atleast"] = "any",
    atleast: int = 3,
    view: View | None = _VIEWS["log2fc"],
) -> Callable[[DataFrame], DataFrame]:

    def filter_func(df: DataFrame) -> DataFrame:
        selected_columns = (
            columns if columns is not None else df.columns.to_list()
        )  # noqa: E501
        passed = df[selected_columns].abs().gt(threshold)

        if how == "any":
            mask = passed.any(axis=1)
        elif how == "all":
            mask = passed.all(axis=1)
        elif how == "atleast":
            mask = passed.sum(axis=1) >= atleast
        else:
            raise ValueError(
                f"Invalid 'how' parameter: {how}. Use 'any', 'all', or 'atleast'."  # noqa: E501
            )

        return df[mask]

    return filter_func


def filter_morphs(filter_specs: Mapping[str, Any]) -> list[Morph]:
    morphs: list[Morph] = []
    for key, spec in filter_specs.items():
        morphs.append(morph(**spec))
    return morphs


def de_filter(
    sample_id: list[str],
    log2_view: View | None = None,
    fdr_view: View | None = None,
    alpha: float = 0.05,
    logthreshold: float = 1,
) -> list[Morph]:
    log2_view = log2_view or View(metric=("Log2FC"), samples=sample_id)
    log_fdr_view = fdr_view or View(metric=("FDR"), samples=sample_id)

    return [
        Morph.from_callable(
            filter_to_threshold(
                threshold=alpha,
                how="any",
                view=log_fdr_view,
                greater_than=False,
            )
        ),
        Morph.from_callable(
            filter_abs_threshold(
                columns="log2fc",
                threshold=logthreshold,
                how="any",
                view=log2_view,
            )
        ),
    ]


def as_table_source(x: Element) -> TableSource:
    if isinstance(x.primary, TableArtifact):
        return x.primary
    return PathTableSource(x.primary.resolve())


def artifacts_for_mode(path: Path, mode: str = "both") -> dict[str, Path]:
    """
    deprecated, we don't usde mode any more, but we keep it for backward
    compatibility.
    """
    if mode == "tsv":
        return {"tsv": path.with_suffix(".tsv")}
    elif mode == "parquet":
        return {"parquet": path.with_suffix(".parquet")}
    elif mode == "both":
        return {
            "tsv": path.with_suffix(".tsv"),
            "parquet": path.with_suffix(".parquet"),
        }
    else:
        raise ValueError(f"Invalid mode: {mode}")


def bundle_morphs(
    schema: ColumnSchema, morphs: list[Morph]
) -> tuple[ColumnSchema, list[dict[str, str]]]:
    running_schema = schema
    select_columns: list[list[str]] = []
    renames: list[dict[str, str]] = []
    drops: list[list[str]] = []

    for m in morphs:
        in_cols = (
            running_schema.select(view=m.view_in)
            if m.view_in
            else list(running_schema)  # noqa: E501
        )
        print("Morph", m)
        if m.view_out is None:
            renames.append({})
            select_columns.append(in_cols)
            continue
        if m.create_new:
            # we create completely new columns based on the view alone
            derived_schema = running_schema.derive_from_view(m.view_out)
            renames.append({})  # no need to rename
            drops.append([]) if m.add else drops.append(in_cols)
        else:
            # we modify the inputz columns create new column names
            derived_schema, rename = running_schema.derive(
                in_cols,
                m.view_out,
                m.append_metric,
            )
        if m.add:
            running_schema = running_schema.merge(derived_schema)
        else:
            running_schema = running_schema.drop(in_cols).merge(derived_schema)
        renames.append(rename)
        drops.append([]) if m.add else drops.append(in_cols)
    return running_schema, renames, drops


def divide(
    numerator: str, denominator: str, *, views: Mapping[str, View]
) -> MultiMorph:

    @morph
    def morph_divide(dataframes: Mapping[str, DataFrame]) -> DataFrame:
        num_df = dataframes[numerator]
        denom_df = dataframes[denominator]
        result_df = num_df.divide(denom_df)
        rename = {}
        for col in num_df.columns:
            old = ColumnTag.decode(col)
            new_tag = ColumnTag(
                metric=old.metric + ("Norm",),
                sample_id=old.sample_id,
                pipeline=old.pipeline,
            )
            rename[col] = new_tag.encode()
        result_df = result_df.rename(columns=rename)
        return result_df

    return morph_divide


################################################################################
# tables
################################################################################


class Tables:
    """
    A class for dataframe operations, including filtering, merging, and
    transforming dataframes for the elements framework.
    """

    def __init__(self):
        self.default_output_spec = OutputSpec(
            ext="parquet", additional_extensions=["tsv"]
        )

    @property
    def default_dir(self) -> Path:
        """
        Default directory for table outputs.
        """
        return Path("results/tables")

    def _finalize(
        self,
        *,
        source: Element,
        tag: ElementTag,
        artifacts: ArtifactSet,
        run_body: Callable[[], None],
        determinants: tuple[str, ...] | None,
        depends_on: tuple[Callable, ...] = (),
    ) -> Element:
        runner = depends(*depends_on)(run_body) if depends_on else run_body
        key, name = generate_element_key_name(tag, "Tables")
        return Element(
            key,
            runner,
            tag=tag,
            determinants=determinants,
            inputs=(source.primary.resolve(),),
            artifacts=artifacts,
            pres=(source,),
            name=name,
        )

    ############################################################################
    # Table Elements
    ############################################################################

    @element
    def transform(
        self,
        source: Element,
        *morphs: Morph,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        # column_schema: ColumnSchema | None = None,
        view: View | None = None,
        index_column: str | None = None,
        output_spec: OutputSpec | None = None,
        params: Params | None = None,
    ) -> Element:
        """
        Transform a dataframe element using the provided morph function. For
        convenience, this is intended only for column transformations, and not
        for row filtering or other operations that would change the number of
        rows in the dataframe. For those operations, use the `filter` method
        instead.

        This method applies the given morph function to the dataframe contained
        within the element, producing a new transformed element.

        Parameters
        ----------
        source : Element
            The input element containing the dataframe to be transformed.
        morph : Callable[[DataFrame], DataFrame]
            A function that takes a DataFrame and returns a transformed
            DataFrame.
        """
        tag = from_prior(
            source.tag,
            tag,
            root=root,
            state=State.TRANSFORMED,
            method=Method.TABLES,
        )

        table_source = as_table_source(source)
        infile = table_source.resolve()
        # running_schema, renames, drops = bundle_morphs(
        #     table_source.schema, list(morphs)
        # )
        output_spec = output_spec or self.default_output_spec
        artifacts, output = ArtifactSet.generate_file_artifacts(
            tag=tag,
            infile=infile,
            spec=output_spec,
            # column_schema=None,  # column_schema or running_schema,
            index_column=index_column or table_source.index_column,
        )  # in another function

        @depends(frame_callable, *(m.fn for m in morphs))
        def __run():
            df = table_source.view(view)
            # columns = table_source.schema.select(view=view) if view else df.columns   # noqa: E501
            # df = df[columns]  # select columns for this view, if any
            for m in morphs:
                # for m, rename in zip(morphs, renames, drops):
                # cols = [c for c in rename]  #
                # result = m.fn(df[cols]).rename(columns=rename)
                df = m.fn(df)  # .rename(columns=rename)
                # df = pd.concat(
                #     [df, result], axis=1
                # )

            print(df.head())
            write_frames(df, output)

        return self._finalize(
            source=source,
            tag=tag,
            artifacts=artifacts,
            run_body=__run,
            determinants=(str(params),) if params else None,
            depends_on=(frame_callable, *(m.fn for m in morphs), bundle_morphs),
        )

    @element
    def filter(
        self,
        source: Element,
        *morphs: Morph,
        root: str | None = None,
        view: View | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        output_spec: OutputSpec | None = None,
        params: Params | None = None,
    ) -> Element:
        tag = from_prior(
            source.tag,
            tag,
            root=root or source.tag.root,
            state=State.FILTER,
            method=Method.TABLES,
        )
        params = params or Params()

        # default morph if none provided
        if not morphs:
            view = _VIEWS["escd"] if view is None else view

            morphs = (
                Morph.from_callable(
                    filter_to_threshold(**params.to_dict()),
                ),
            )

        infile = source.primary.resolve()
        artifacts, output = ArtifactSet.generate_file_artifacts(
            tag=tag,
            infile=infile,
            spec=output_spec or self.default_output_spec,
            # column_schema=source.primary.column_schema,
            index_column=source.primary.index_column,
        )

        @depends(frame_callable, *(m.fn for m in morphs), bundle_morphs)
        def __run():
            df = source.primary.view(view)
            for morph in morphs:
                df = morph.fn(df)
            write_frames(df, output)

        return self._finalize(
            source=source,
            tag=tag,
            artifacts=artifacts,
            run_body=__run,
            determinants=(str(params),) if params else None,
            depends_on=(frame_callable, *(m.fn for m in morphs), bundle_morphs),
        )

    @element
    def combine(
        self,
        source: Element,
        others: Element | Iterable[Element],
        *morphs: MultiMorph,
        root: str | None = None,
        index_column: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        views: Mapping[str, View] | None = None,
        output_spec: OutputSpec | None = None,
        params: Params | None = None,
    ) -> Element:
        """
        Combine multiple DataFrame elements into a single DataFrame element
        using the provided morph function.

        This method applies the given morph function to the dataframe contained
        within the element, producing a new transformed element. This should
        also be usable for merging, concatenation and joining.

        Parameters
        ----------
        source : Element
            The input element containing the dataframe to be transformed.
        others : Element | Iterable[Element]
            Other elements to be combined with the source element.
        morphs : MultiMorph
            Multiple morph functions that take a mapping of DataFrames and
            return a transformed mapping of DataFrames.
        root : str | None, optional
            Root name for the combined element, by default None
        index_column : str | None, optional
            Index column for the combined element, by default None
        tag : PartialElementTag | ElementTag | None, optional
            Tag for the combined element, by default None
        views : Mapping[str, View] | None, optional
            Views to select columns for each of the input DataFrames, by default
            None (select all columns).
        output_spec : OutputSpec | None, optional
            Output specification for the combined element, by default None
        params : Params | None, optional
            Parameters for the morph functions, by default None
        """
        views = views or {}
        tag = from_prior(
            source.tag,
            tag,
            root=root,
            state=State.COMBINED,
            method=Method.TABLES,
        )
        sources = {
            source.root: as_table_source(source),
            **{
                s.root: as_table_source(s)
                for s in others
                if isinstance(others, Iterable)
            },
        }
        infiles = {k: s.resolve() for k, s in sources.items()}
        output_spec = output_spec or self.default_output_spec
        artifacts, output = ArtifactSet.generate_file_artifacts(
            tag=tag,
            infile=infiles,
            spec=output_spec,
            # column_schema=None,  # column_schema or running_schema,
            index_column=index_column or sources[source.root].index_column,
        )  # in another function

        @depends(frame_callable, *(m.fn for m in morphs))
        def __run():
            table_sources = {
                k: s.view(views.get(k)) for k, s in sources.items()
            }  # noqa: E501
            inputs = table_sources
            for m in morphs:
                inputs = m.fn(inputs)

            for (result_key,) in inputs.items():
                if result_key not in artifacts:
                    raise ValueError(
                        f"Result key '{result_key}' not found in artifacts."
                    )
                output = artifacts[result_key]
                df = inputs[result_key]
                write_frames(df, output)

        key, name = generate_element_key_name(
            tag, "Tables", subroutine="combine"
        )  # noqa: E501
        return Element(
            key,
            __run,
            tag=tag,
            determinants=tuple(m.signature() for m in morphs),
            inputs=tuple(infiles.values()),
            artifacts=artifacts,
            pres=(
                (source, *others.values())
                if isinstance(others, dict)
                else (source, *others)
            ),
            name=name,
        )

    ############################################################################
    # From genes
    ############################################################################

    def fromgenes(
        self,
        source: FileElement,
        sample_names: list[str],
        *,
        # column_schema: ColumnSchema | None = None,
        index_column: str | None = "gene_stable_id [META]",
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        output_spec: OutputSpec | None = None,
    ):
        schema, rename = ColumnSchema.from_legacy(
            read_schema(source.primary), sample_names
        )
        # column_schema = column_schema or schema
        gene_cols = [
            "gene_stable_id",
            "name",
            "chr",
            "start",
            "stop",
            "strand",
            "tss",
            "tes",
            "biotype",
            "transcript_stable_ids",
        ]

        def __rename_legacy(df: DataFrame) -> DataFrame:
            for col in gene_cols:
                df[col] = df[col].astype("str")
            if rename:
                df = df.rename(columns=rename)
                for col in df.columns:
                    print("col", col)
            if index_column:
                df = df.set_index(index_column)
            return df

        morph = Morph(fn=__rename_legacy)

        return self.transform(
            source,
            morph,
            root=root,
            tag=tag,
            # column_schema=column_schema,
            index_column=index_column,
            output_spec=output_spec,
        )

    @element
    def fromfile(
        self,
        path: Path,
        *,
        index_column: str | None = None,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
    ) -> Element:
        tag = PartialElementTag(
            root=root or Path(path).stem,
            state=State.INCOMING,
            method=Method.TABLES,
        ).merge(tag)
        # column_schema = ColumnSchema.build(read_schema(path))
        source = TableArtifact(
            path=Path(path),
            index_column=index_column,
        )  # column_schema=ecolumn_schema
        return FileElement(
            source,
            tag=tag,
            root=root,
        )


################################################################################
# convenience functions
################################################################################


def frame_callable(
    primary: TableArtifact,
    *,
    view: View | None = None,
    additional_filter: Callable[[DataFrame], DataFrame] | None = None,
) -> Callable[[], tuple[DataFrame, DataFrame]]:

    def __load() -> tuple[DataFrame, DataFrame]:
        cols = primary.column_schema.select(view)
        if not cols:
            raise ValueError(f"No columns match view={view!r}")
        df = primary.frame  # via read_frame
        if additional_filter is not None:
            df = additional_filter(df)

        return df, df[cols]

    return __load


def drop_retain_rename(
    select: list[str] | None = None,
    retains: bool = True,
    rename: Mapping[str, str] | None = None,
) -> Callable[[pd.DataFrame], pd.DataFrame]:
    """
    Returns a function that selects or drops columns based on the provided list.
    If `retains` is True, only the columns in the `select` list are retained.
    If `retains` is False, the columns in the `select` list are dropped.
    If a rename mapping is provided, it will also rename the columns according
    to the mapping.
    """

    def _drop_frqs(df: pd.DataFrame) -> pd.DataFrame:
        if rename:
            df = df.rename(columns=rename)
        if select:
            selected_columns = [
                col
                for col in df.columns.tolist()
                if any(sub in col for sub in select)  # noqa: E501
            ]
            if retains:
                assert all(
                    col in df.columns for col in selected_columns
                ), "Some columns in 'select' are not in the DataFrame"
                df = df[[col for col in df.columns if col in selected_columns]]
            else:
                assert all(
                    col in df.columns for col in selected_columns
                ), "Some columns in 'select' are not in the DataFrame"
                df = df[
                    [col for col in df.columns if col not in selected_columns]
                ]  # noqa: E501
        return df

    return _drop_frqs


def validate_merge_columns(
    left_df: pd.DataFrame,
    right_df: pd.DataFrame,
    on_left: list[str],
    on_right: list[str],
) -> None:
    left_cols = set(left_df.columns)
    right_cols = set(right_df.columns)
    join_cols = set(on_right + on_left)
    conflicts = (left_cols & right_cols) - join_cols
    if conflicts:
        raise ValueError(
            f"Column name conflicts detected during merge: {conflicts}. "
            f"Consider renaming these columns in one of the DataFrames before merging."  # noqa: E501
        )


def append_label(label: str) -> Callable[[str], str]:
    def _transform_name(name: str) -> str:
        if "(" not in name or not name.endswith(")"):
            return f"{name}.{label}"

        base, suffix = name.rsplit("(", 1)
        return f"{base.rstrip()}.{label} ({suffix})"

    return _transform_name


@morph(view_out=View(metric=("Z")), add=True)
def zscore(df: pd.DataFrame, rename, drop) -> pd.DataFrame:
    """
    Computes the z-score for specified columns in a DataFrame.
    The z-score is calculated as (x - mean) / std for each value x in column.
    """
    columns = df.columns
    rename = {}
    for col in columns:
        old = ColumnTag.decode(col)
        new_tag = ColumnTag(
            metric=old.metric + ("Z",),
            sample_id=old.sample_id,
            pipeline=old.pipeline,  # noqa: E501
        )
        rename[col] = new_tag.encode()
    df = df.rename(columns=rename)
    for col in df.columns:
        mean = df[col].mean()
        std = df[col].std()
        df[col] = (df[col] - mean) / std
    return df
