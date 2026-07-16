from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property, wraps
from pathlib import Path
from typing import Any, Callable, Iterable, Literal, Mapping, Protocol

import pandas as pd  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

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
from mmalignments.services.io import (
    read_frame,
    read_schema,
    write_frame,
    write_frames,
)

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
        cols = (
            list(self.schema.select(view)) if view else list(self.frame.columns)
        )  # noqa: E501
        return self.frame[cols]  # type: ignore[return-value]


################################################################################
# Morphs
################################################################################


class MorphProtocol(Protocol):
    fn: Callable
    view_in: View | None
    params: Mapping[str, str | float | int | bool]


@dataclass(frozen=True)
class Morph:
    fn: Callable[
        [
            DataFrame,
        ],
        DataFrame,
    ]
    view_in: View | None = None
    params: Mapping[str, str | float | int | bool] = field(default_factory=dict)
    # view_out: View | None = None
    # create_new: bool = False
    # add: bool = False
    # append_metric: bool = True

    @classmethod
    def from_callable(
        cls,
        fn: Callable,
        *,
        view: View | None = None,
        params: dict[str, Any] | None = None,
    ) -> Morph:  # noqa: E501
        return cls(fn=fn, view_in=view, params=params or {})

    def __str__(self):
        return f"Morph(fn={self.fn.__name__}, view_in={self.view_in}, params={self.params})"  # noqa: E501

    def signature(self):
        params_hash = stable_hash(self.params) + stable_hash(self.view_in)
        return function_hash(self.fn) + params_hash


def morph(
    _fn: Callable | None = None,
    *,
    view_in: View | None = None,
    params: Mapping[str, str | float | int | bool] | None = None,
):
    """Decorator turning a function into a Morph."""

    def decorator(fn: Callable[..., DataFrame]) -> Morph:
        @wraps(fn)
        def wrapped(*args, **kwargs):
            return fn(*args, **kwargs)

        return Morph(
            fn=wrapped,
            view_in=view_in,
            params=params or {},
        )

    if _fn is None:
        return decorator

    return decorator(_fn)


@dataclass(frozen=True)
class MultiMorph:
    fn: Callable[[Mapping[str, DataFrame]], Mapping[str, DataFrame]]
    views: Mapping[str, View] | None = None
    params: Mapping[str, Mapping[str, Any]] | None = field(default_factory=dict)

    def __str__(self):
        return f"MultiMorph(fn={self.fn.__name__}, views={self.views})"  # noqa: E501

    def signature(self):
        view_hash = stable_hash(self.views) if self.views else ""
        params_hash = stable_hash(self.params) if self.params else ""
        return function_hash(self.fn) + params_hash + view_hash

    @classmethod
    def from_callable(
        cls,
        fn: Callable[[Mapping[str, DataFrame]], Mapping[str, DataFrame]],  # noqa: E501
        *,
        views: Mapping[str, View] | None = None,
    ) -> MultiMorph:

        morph_views = getattr(fn, "morph_views", views)
        return cls(fn=fn, views=morph_views if morph_views else None)  # noqa: E501


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
) -> Morph:

    def filter_func(df: DataFrame) -> DataFrame:
        count_columns = columns if columns is not None else df.columns.to_list()
        passed = (
            df[count_columns].gt(threshold)
            if greater_than
            else df[count_columns].lt(threshold)
        )

        if how == "any":
            mask = passed.any(axis=1)  # type: ignore[return-value]
        elif how == "all":
            mask = passed.all(axis=1)  # type: ignore[return-value]
        elif how == "atleast":
            mask = passed.sum(axis=1) >= atleast  # type: ignore[return-value]
        else:
            raise ValueError(
                f"Invalid 'how' parameter: {how}. Use 'any', 'all', or 'atleast'."  # noqa: E501
            )

        return df.loc[mask, :]  # type: ignore[return-value]

    return Morph.from_callable(
        filter_func,
        view=view,
        params={
            "threshold": threshold,
            "how": how,
            "atleast": atleast,
            "columns": columns,
            "greater_than": greater_than,
        },
    )


def filter_abs_threshold(
    columns: list[str] | None = None,
    threshold: int | float = 1,
    how: Literal["any", "all", "atleast"] = "any",
    atleast: int = 3,
    view: View | None = View(metric="log2FC", sample_id="Miss_vs_Non"),
) -> Morph:

    def filter_func(df: DataFrame) -> DataFrame:
        selected_columns = columns or df.columns.to_list()
        if view:
            selected_columns = ColumnTag.select_from_view(
                selected_columns, view
            )  # noqa: E501

        passed = df[selected_columns].abs().gt(threshold)

        if how == "any":
            mask = passed.any(axis=1)  # type: ignore[return-value]
        elif how == "all":
            mask = passed.all(axis=1)  # type: ignore[return-value]
        elif how == "atleast":
            mask = passed.sum(axis=1) >= atleast  # type: ignore[return-value]
        else:
            raise ValueError(
                f"Invalid 'how' parameter: {how}. Use 'any', 'all', or 'atleast'."  # noqa: E501
            )

        return df[mask]  # type: ignore[return-value]

    return Morph.from_callable(
        filter_func,
        params={
            "columns": columns,
            "threshold": threshold,
            "how": how,
            "atleast": atleast,
            "view": view,
        },
    )


def filter_morphs(filter_specs: Mapping[str, Any]) -> list[Morph]:
    morphs: list[Morph] = []
    for key, spec in filter_specs.items():
        morphs.append(Morph(**spec))
    return morphs


def de_filter(
    columns_log: tuple[str] | None = None,
    columns_fdr: tuple[str] | None = None,
    log2_view: View | None = None,
    fdr_view: View | None = None,
    alpha: float = 0.05,
    logthreshold: float = 1,
) -> list[Morph]:
    log2_view = log2_view or View(metric=("Log2FC"))
    log_fdr_view = fdr_view or View(metric=("FDR"))
    return [
        filter_to_threshold(
            threshold=alpha,
            how="any",
            view=log_fdr_view,
            greater_than=False,
        ),
        filter_abs_threshold(
            columns=None,
            threshold=logthreshold,
            how="any",
            view=log2_view,
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


# def bundle_morphs(
#     schema: ColumnSchema, morphs: list[Morph]
# ) -> tuple[ColumnSchema, list[dict[str, str]]]:
#     running_schema = schema
#     select_columns: list[list[str]] = []
#     renames: list[dict[str, str]] = []
#     drops: list[list[str]] = []

#     for m in morphs:
#         in_cols = (
#             running_schema.select(view=m.view_in)
#             if m.view_in
#             else list(running_schema)  # noqa: E501
#         )
#         print("Morph", m)
#         if m.view_out is None:
#             renames.append({})
#             select_columns.append(in_cols)
#             continue
#         if m.create_new:
#             # we create completely new columns based on the view alone
#             derived_schema = running_schema.derive_from_view(m.view_out)
#             renames.append({})  # no need to rename
#             drops.append([]) if m.add else drops.append(in_cols)
#         else:
#             # we modify the inputz columns create new column names
#             derived_schema, rename = running_schema.derive(
#                 in_cols,
#                 m.view_out,
#                 m.append_metric,
#             )
#         if m.add:
#             running_schema = running_schema.merge(derived_schema)
#         else:
#            running_schema = running_schema.drop(in_cols).merge(derived_schema)
#         renames.append(rename)
#         drops.append([]) if m.add else drops.append(in_cols)
#     return running_schema, renames, drops


def divide(
    numerator: str,
    denominator: str,
    denominator_column: str,
    views: Mapping[str, View] | None = None,
    *,
    rename_fn: Callable[[list[str]], dict[str, str]] | None = None,
) -> MultiMorph:

    def morph_divide(
        dataframes: Mapping[str, DataFrame],
    ) -> Mapping[str, DataFrame]:  # noqa: E501
        num_df = dataframes[numerator]
        denom_df = dataframes[denominator]

        result_df = num_df.divide(
            denom_df[denominator_column].reindex(num_df.columns), axis="columns"
        )
        if rename_fn:
            rename = rename_fn(result_df.columns.to_list())
            result_df = result_df.rename(columns=rename)
        rename = {}
        for col in num_df.columns:
            old = ColumnTag.decode(col)
            new_metric = (
                old.metric + ("Norm",)
                if isinstance(old.metric, tuple)
                else (old.metric, "Norm")
            )

            new_tag = ColumnTag(
                sample_id=old.sample_id,
                pipeline=old.pipeline,
                metric=new_metric,
            )
            rename[col] = new_tag.encode()
        result_df = result_df.rename(columns=rename)
        return {"primary": result_df}

    return MultiMorph.from_callable(morph_divide, views=views)


################################################################################
# tables
################################################################################


class Tables:
    """
    A class for dataframe operations, including filtering, merging, and
    transforming dataframes for the elements framework.
    """

    def __init__(self):
        pass


    @property
    def default_dir(self) -> Path:
        """
        Default directory for table outputs.
        """
        return Path("results/tables")

    @cached_property
    def default_output_spec(self) -> OutputSpec:
        ret = OutputSpec(
            outdir=self.default_dir,
            ext="parquet",
            additional_extensions=("tsv",),
        )
        return ret
    # def _finalize(
    #     self,
    #     *,
    #     source: Element,
    #     tag: ElementTag,
    #     artifacts: ArtifactSet,
    #     run_body: Callable[[], None],
    #     determinants: tuple[str, ...] | None,
    #     depends_on: tuple[Callable, ...] = (),
    # ) -> Element:
    #     # runner = depends(*depends_on)(run_body) if depends_on else run_body
    #     key, name = generate_element_key_name(tag, "Tables")
    #     return Element(
    #         key,
    #         runner,
    #         tag=tag,
    #         determinants=determinants,
    #         inputs=(source.primary.resolve(),),
    #         artifacts=artifacts,
    #         pres=(source,),
    #         name=name,
    #     )

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
            default_dir=infile.parent,
            spec=output_spec,
            # column_schema=None,  # column_schema or running_schema,
            index_column=index_column or table_source.index_column,
        )  # in another function

        @depends(frame_callable, *(m.fn for m in morphs))
        def __run():
            df = table_source.view(view)
            for m in morphs:
                df = m.fn(df)  # .rename(columns=rename)

            write_frames(df, output)

        key, name = generate_element_key_name(
            tag, "Tables", subroutine="transform"
        )  # noqa: E501
        return Element(
            key,
            __run,
            tag=tag,
            determinants=(str(params),) if params else None,
            inputs=(source.primary.resolve(),),
            artifacts=artifacts,
            pres=(source,),
            name=name,
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
            morphs = (filter_to_threshold(**params.to_dict()),)

        infile = source.primary.resolve()
        artifacts, output = ArtifactSet.generate_file_artifacts(
            tag=tag,
            default_dir=infile.parent,
            spec=output_spec or self.default_output_spec,
            # column_schema=source.primary.column_schema,
            index_column=source.primary.index_column,
        )

        # @depends(frame_callable, *(m.fn for m in morphs), bundle_morphs)
        def __run():
            filter_df = source.primary.view(view)
            df = source.primary.view()
            for morph in morphs:
                index = morph.fn(filter_df).index
                df = df.loc[index, :]
            write_frames(df, output)

        # return self._finalize(
        #     source=source,
        #     tag=tag,
        #     artifacts=artifacts,
        #     run_body=__run,
        #     determinants=(str(params),) if params else None,
        #   depends_on=(frame_callable, *(m.fn for m in morphs), bundle_morphs),
        # )
        key, name = generate_element_key_name(
            tag, "Tables", subroutine="filter"
        )  # noqa: E501
        return Element(
            key,
            __run,
            tag=tag,
            determinants=(str(params),) if params else None,
            inputs=(source.primary.resolve(),),
            artifacts=artifacts,
            pres=(source,),
            name=name,
        )

    @element
    def combine(
        self,
        sources: Mapping[str, Element],
        *morphs: MultiMorph,
        artifact_keys: Mapping[str, str] | None = None,
        first_key: str | None = None,
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
        sources : Mapping[str, Element]
            The input elements containing the dataframes to be transformed and
            combined.
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
        artifact_keys = artifact_keys or {}
        first_element = (
            sources[first_key] if first_key else next(iter(sources.values()))
        )
        views = views or {}
        tag = from_prior(
            first_element.tag,
            tag,
            root=root,
            state=State.COMBINED,
            method=Method.TABLES,
        )
        input_paths = {
            key: s.artifacts[artifact_keys.get(key, "primary")].resolve()
            for key, s in sources.items()
        }
        infile = first_element.primary.resolve()
        output_spec = output_spec or self.default_output_spec
        artifacts, _ = ArtifactSet.generate_file_artifacts(
            tag=tag,
            default_dir=infile.parent,
            spec=output_spec,
            index_column=index_column or first_element.primary.index_column,
        )  # in another function

        @depends(frame_callable)
        def __run():
            inputs = {}
            for key, source in sources.items():
                inputs[key] = source.artifacts[
                    artifact_keys.get(key, "primary")
                ].view(  # noqa: E501
                    views.get(key)
                )
            # apply the morphs
            for m in morphs:
                inputs = m.fn(inputs)

            results = (
                {artifacts.primary_name: inputs}
                if isinstance(inputs, DataFrame)
                else inputs
            )
            for result_key, df in results.items():
                if result_key not in artifacts:
                    raise ValueError(
                        f"Result key '{result_key}' not found in artifacts."
                    )
                out = artifacts.primary.resolve()
                write_frame(df, out)

        key, name = generate_element_key_name(
            tag, "Tables", subroutine="combine"
        )  # noqa: E501
        return Element(
            key,
            __run,
            tag=tag,
            determinants=tuple(m.signature() for m in morphs),
            inputs=tuple(s.resolve() for s in input_paths.values()),
            artifacts=artifacts,
            pres=tuple(sources.values()),
            name=name,
        )

    @element
    def join(
        self,
        *sources: tuple[Element, str] | Element,
        on: str | None = None,
        how: Literal[
            "left",
            "right",
            "outer",
            "inner",
            "cross",
            "left_anti",
            "right_anti",  # noqa: E501
        ] = "left",
        views: Mapping[str, View | None] | None = None,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        output_spec: OutputSpec | None = None,
    ) -> Element:
        """
        Combine multiple DataFrame elements into a single DataFrame element
        using the provided morph function.

        This method applies the given morph function to the dataframe contained
        within the element, producing a new transformed element. This should
        also be usable for merging, concatenation and joining.

        Parameters
        ----------
        sources : Element
            an arbitrary number of Element instances.
        on : str, optional
            Column name to join on, by default "__index"
        how : str, optional
            Type of join to perform, by default "inner"
        root : str | None, optional
            Root name for the combined element, by default None
        tag : PartialElementTag | ElementTag | None, optional
            Tag for the combined element, by default None
        output_spec : OutputSpec | None, optional
            Output specification for the combined element, by default None

        Returns
        -------
        Element
            A new Element instance representing the joined DataFrame.
        """

        def get_element_and_key(
            source: tuple[Element, str] | Element,
        ) -> tuple[Element, str]:
            if isinstance(source, tuple):
                return source
            else:
                return source, "primary"

        views = views or {}
        output_spec = output_spec or self.default_output_spec
        first_element, fkey = get_element_and_key(sources[0])
        fkey = fkey or "primary"
        tag = from_prior(
            first_element.tag,
            tag,
            root=root or first_element.tag.root,
            state=State.COMBINED,
            method=Method.TABLES,
        )
        elements = []
        names_in_order = []
        pre = []
        for source in sources:
            el, key = get_element_and_key(source)
            key = key or "primary"
            pre.append(el)
            names_in_order.append(el.name)
            elements.append(el)

        artifacts, output = ArtifactSet.generate_file_artifacts(
            tag=tag,
            default_dir=first_element.artifacts[fkey].resolve().parent,
            spec=output_spec,
            index_column=first_element.artifacts[fkey].index_column,
        )  # in another function

        def __load_frame(source: Path | TableArtifact) -> DataFrame:
            if isinstance(source, TableArtifact):
                return source.view()
            else:
                return read_frame(source)

        @depends(__load_frame)
        def __run():
            df: DataFrame | None = None
            for ii, key in enumerate(names_in_order):
                el = elements[ii]
                view = views.get(key)
                other_df: DataFrame = el.artifacts["primary"].view(view)
                if df is None:
                    df = other_df
                else:
                    df = df.join(other_df, on=on, how=how)
            if df is None:
                raise ValueError("Cannot join empty input artifacts.")
            write_frames(df, output)

        key, name = generate_element_key_name(
            tag, "Tables", subroutine="join"
        )  # noqa: E501
        return Element(
            key,
            __run,
            tag=tag,
            determinants=(str(on), str(how)),
            artifacts=artifacts,
            pres=tuple(pre),
            name=name,
        )

    def long(
        self,
        source: FileElement,
        *,
        sample_filter: Callable[[str | None], bool] | None = None,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        index_column: str | None = "gene_stable_id [META]",
        output_spec: OutputSpec | None = None,
    ):
        """
        Convert a wide-format table to long-format.

        Parameters
        ----------
        source : FileElement
            The input FileElement containing the wide-format table.
        sample_filter : Callable[[ColumnTag], bool] | None, optional
            A function to filter sample columns, by default None
        root : str | None, optional
            The root tag for the output element, by default None
        tag : PartialElementTag | ElementTag | None, optional
            The tag for the output element, by default None
        index_column : str | None, optional
            The column to use as the index, by default "gene_stable_id [META]"
        output_spec : OutputSpec | None, optional
            The output specification, by default None

        Returns
        -------
        Element
            The output Element containing the long-format table.
        """
        tag = from_prior(
            source.tag,
            tag,
            root=root,
            state=State.TRANSFORMED,
            method=Method.TABLES,
            param="long",
        )
        return self.transform(
            source,
            self.to_long(
                index=index_column, sample_filter=sample_filter
            ),  # pivot table
            root=root,
            tag=tag,
            index_column=index_column,
            output_spec=output_spec,
        )

    def to_long(
        self,
        index: str | None = "gene_stable_id [META]",
        sample_filter: Callable[[str | None], bool] | None = None,
    ) -> Morph:
        def sample_filter_default(sample_id: str | None) -> bool:
            return sample_id is not None

        def __transform(df: DataFrame) -> DataFrame:
            if index in df.columns:
                df = df.set_index(index)

            tags = {c: ColumnTag.decode(c) for c in df.columns}
            is_sample = sample_filter or sample_filter_default
            for c, tag in tags.items():
                print(tag.sample_id)
            # raise NotImplementedError("sample_filter is not implemented yet.")
            sample_cols = [
                c for c, tag in tags.items() if is_sample(tag.sample_id)
            ]  # columns with samples in ()
            meta_cols = [
                c for c, tag in tags.items() if tag.sample_id is None
            ]  # columns with META pipeline

            # vectorized lookup maps statt apply(pd.Series)
            sample_id_map = {c: tags[c].sample_id for c in sample_cols}
            pipeline_map = {c: tags[c].pipeline for c in sample_cols}
            metric_map = {
                c: (
                    ".".join(tags[c].metric)
                    if isinstance(tags[c].metric, tuple)
                    else tags[c].metric
                )
                for c in sample_cols
            }

            long = (
                df[sample_cols]
                .reset_index()
                .melt(
                    id_vars=[df.index.name],
                    var_name="column",
                    value_name="value",  # noqa: E501
                )  # long format
            )
            long["sample_id"] = long["column"].map(sample_id_map)
            long["pipeline"] = long["column"].map(pipeline_map)
            long["metric"] = long["column"].map(metric_map)
            long = long.drop(columns="column")
            # pivot statt pivot_table -> kein Aggregations-Overhead
            long = long.pivot(
                index=[df.index.name, "sample_id", "pipeline"],
                columns="metric",
                values="value",
            ).reset_index()
            if meta_cols:
                long = long.merge(
                    df[meta_cols].reset_index(), on=df.index.name, how="left"
                )
            return long

        return Morph.from_callable(__transform)

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
            mask = df["gene_stable_id"].str.startswith("ENS")
            df = df.loc[mask, :]  # type: ignore[return-value]
            for col in gene_cols:
                df[col] = df[col].astype("str")
            if rename:
                df = df.rename(columns=rename)
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
        cols = primary.schema.select(view)
        if not cols:
            raise ValueError(f"No columns match view={view!r}")
        df = primary.frame  # via read_frame
        if additional_filter is not None:
            df = additional_filter(df)

        return df, df[cols]  # type: ignore[return-value]

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
                df = df[[col for col in df.columns if col in selected_columns]]  # type: ignore[return-value]
            else:
                assert all(
                    col in df.columns for col in selected_columns
                ), "Some columns in 'select' are not in the DataFrame"
                df = df[
                    [col for col in df.columns if col not in selected_columns]  # type: ignore[return-value]
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


@morph
def zscore(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes the z-score for specified columns in a DataFrame.
    The z-score is calculated as (x - mean) / std for each value x in column.
    """
    rename = ColumnTag.rename_columns(
        df.columns.to_list(), view=_VIEWS["zscore"], append=True
    )
    df = df.rename(columns=rename)
    for col in df.columns:
        mean = df[col].mean()
        std = df[col].std()
        df[col] = (df[col] - mean) / std
    return df


@morph
def cpm(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes counts per million (CPM) for specified columns in a DataFrame.
    CPM is calculated as (x / column_sum) * 1e6 for each value x in column.
    """
    rename = ColumnTag.rename_columns(
        df.columns.to_list(), view=_VIEWS["cpm"], append=True
    )
    df = df.rename(columns=rename)
    for col in df.columns:
        column_sum = df[col].sum()
        df[col] = (df[col] / column_sum) * 1e6
    return df


def append_value(columns_values: dict[str, str]) -> Morph:
    """
    Returns a Morph that appends new columns with specified values to a DataFrame.

    Parameters
    ----------
    columns_values : dict[str, str]
        Dictionary mapping column names to the values to be appended.

    Returns
    -------
    Morph
        Morph that appends new columns with specified values to a DataFrame.
    """

    def __append(df: pd.DataFrame) -> pd.DataFrame:
        """
        Appends new columns with specified values to the DataFrame.
        """
        for key, value in columns_values.items():
            if key in df.columns:
                raise ValueError(f"Column '{key}' already exists in the DataFrame.")
            df[key] = value
        return df

    return Morph.from_callable(__append, params=columns_values)


def left_join(
    other: Element,
    *,
    on: str | None = None,
    view: View | None = None,
    artifact: str = "primary",
    how: Literal[
        "left",
        "right",
        "outer",
        "inner",
        "cross",
        "left_anti",
        "right_anti",
    ] = "left",
) -> Morph:
    def __left_join(df: DataFrame) -> DataFrame:
        print(other.artifacts[artifact].view())
        other_df = other.artifacts[artifact].view(view)
        return df.join(other_df, on=on, how=how)

    return Morph.from_callable(
        __left_join, params={"on": on, "view": view, "artifact": artifact, "how": how}
    )


def to_long(
    index: str | None = "gene_stable_id [META]",
    sample_filter: Callable[[str | None], bool] | None = None,
) -> Morph:
    def sample_filter_default(sample_id: str | None) -> bool:
        return sample_id is not None

    def __transform(df: DataFrame) -> DataFrame:
        if index in df.columns:
            df = df.set_index(index)

        tags = {c: ColumnTag.decode(c) for c in df.columns}
        is_sample = sample_filter or sample_filter_default
        for c, tag in tags.items():
            print(tag.sample_id)
        # raise NotImplementedError("sample_filter is not implemented yet.")
        sample_cols = [
            c for c, tag in tags.items() if is_sample(tag.sample_id)
        ]  # columns with samples in ()
        meta_cols = [
            c for c, tag in tags.items() if tag.sample_id is None
        ]  # columns with META pipeline

        # vectorized lookup maps statt apply(pd.Series)
        sample_id_map = {c: tags[c].sample_id for c in sample_cols}
        pipeline_map = {c: tags[c].pipeline for c in sample_cols}
        metric_map = {
            c: (
                ".".join(tags[c].metric)
                if isinstance(tags[c].metric, tuple)
                else tags[c].metric
            )
            for c in sample_cols
        }

        long = (
            df[sample_cols]
            .reset_index()
            .melt(
                id_vars=[df.index.name],
                var_name="column",
                value_name="value",  # noqa: E501
            )  # long format
        )
        long["sample_id"] = long["column"].map(sample_id_map)
        long["pipeline"] = long["column"].map(pipeline_map)
        long["metric"] = long["column"].map(metric_map)
        long = long.drop(columns="column")
        # pivot statt pivot_table -> kein Aggregations-Overhead
        long = long.pivot(
            index=[df.index.name, "sample_id", "pipeline"],
            columns="metric",
            values="value",
        ).reset_index()
        if meta_cols:
            long = long.merge(df[meta_cols].reset_index(), on=df.index.name, how="left")
        return long

    return Morph.from_callable(__transform)


def morph_rank_file(id_column: str = "id") -> Morph:
    def __transform(data: DataFrame) -> DataFrame:
        data[id_column] = data[id_column].str.capitalize()
        return data

    return Morph.from_callable(__transform)


def filter_to_group(
    group: str, values: list[str], *, view: View | None = None
) -> Morph:
    def __filter(df: DataFrame) -> DataFrame:
        if view is not None:
            selected_columns= ColumnTag.select_from_view(df.columns.to_list(), view)
            if selected_columns:
                df = df[selected_columns]   # type: ignore[assignment]
        return df[df[group].isin(values)]   # type: ignore[return-value]

    return Morph.from_callable(__filter, view=view, params={"group": group})
