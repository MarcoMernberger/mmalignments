from abc import ABC, abstractmethod
from dataclasses import dataclass
import itertools
from pathlib import Path
from typing import Any, Callable, Literal, Mapping

import pandas as pd
from pandas import DataFrame, Series
from mmalignments.services.dependencies import depends
from mmalignments.models.elements import (
    CallSpec,
    Element,
    Runnable,
    TableElement,
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
from mmalignments.services.io import read_frame, write_frames


class Morph(ABC):

    @abstractmethod
    def apply(self, df: DataFrame) -> DataFrame:
        pass


@dataclass(frozen=True)
class DropColumns(Morph):
    columns: tuple[str, ...]

    def apply(self, df: DataFrame) -> DataFrame:
        return df.drop(columns=list(self.columns))


@dataclass(frozen=True)
class FilterRows(Morph):
    func: Callable[[DataFrame], DataFrame]

    def apply(self, df: DataFrame) -> DataFrame:
        return self.func(df)


@dataclass(frozen=True)
class AddColumn(Morph):
    name: str
    func: Callable[[DataFrame], Series]

    def apply(self, df: DataFrame) -> DataFrame:
        df = df.copy()
        df[self.name] = self.func(df)
        return df


@dataclass(frozen=True)
class Merge(Morph):
    other: "TableElement"
    kwargs: dict

    def apply(self, df: DataFrame) -> DataFrame:
        other_df = self.other.read()
        return df.merge(other_df, **self.kwargs)


@dataclass(frozen=True)
class TableSeed:

    source: "TableElement"
    morphs: tuple[Morph, ...] = ()

    def _append(self, morph: Morph) -> "TableSeed":
        return TableSeed(
            source=self.source,
            morphs=self.morphs + (morph,),
        )

    def drop(self, *columns: str) -> "TableSeed":
        return self._append(DropColumns(tuple(columns)))

    def filter(self, func) -> "TableSeed":
        return self._append(FilterRows(func))

    def add_column(self, name: str, func) -> "TableSeed":
        return self._append(AddColumn(name, func))

    def merge(self, other: "TableElement", **kwargs) -> "TableSeed":
        return self._append(Merge(other, kwargs))

    def evaluate(self) -> pd.DataFrame:

        df = self.source.read()

        for morph in self.morphs:
            df = morph.apply(df)

        return df

    @element
    def materialize(
        self,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        mode: str = "both",
    ) -> "TableElement":
        tag = from_prior(
            self.source.tag,
            tag,
            state=State.TRANSFORMED,
            method=Method.TABLES,
        )
        key, name = generate_element_key_name(tag, "Tables", "materialize")
        out_dir = Path(outdir or self.source.file.parent)
        filename = filename or tag.default_output
        outfile = out_dir / filename

        def run():
            df = self.evaluate()
            write_frames(df, outfile, mode=mode)

        callspec = CallSpec(path=("tables", "materialize")).render()
        runner = Runnable(run, display=callspec)

        return TableElement(
            key=key,
            run=runner,
            tag=tag,
            pres=(self.source,),
            artifacts=None,
            name=name,
            index=None,
        )


class Tables:

    @classmethod
    def artifacts_for_mode(cls, path: Path, mode: str = "both") -> dict[str, Path]:
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

    @element
    def concat(
        self,
        *elements: Element,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        mode: Literal["tsv", "parquet", "both"] = "both",
        params: Params | None = None,
    ) -> Element:

        first_element = elements[0]
        root = root or first_element.tag.root
        tag = from_prior(
            first_element.tag,
            state=State.MERGED,
            method=Method.TABLES,
        )
        outdir = Path(outdir or first_element.file.parent)
        filename = filename or tag.default_output
        output_file = Path(outdir) / filename
        frame_paths = {element.tag.root: element.file for element in elements}
        artifacts = Tables.artifacts_for_mode(output_file, mode=mode)
        run = self.concat_tables(frame_paths, output_file, mode=mode, params=params)

        key, name = generate_element_key_name(tag, "Tables")

        return Element(
            key,
            run,
            tag=tag,
            determinants=params.determinants() if params else None,
            inputs=tuple(frame_paths.values()),
            artifacts=artifacts,
            pres=tuple(elements),
            name=name,
        )

    def concat_tables(
        self,
        paths: dict[str, Path],
        output_path: Path,
        *,
        mode: Literal["tsv", "parquet", "both"] = "both",
        fill_columns: dict[str, Any] | None = None,
        params: Params | None = None,
    ) -> Runnable:
        params = params or Params()

        def __run():
            # Concatenate the DataFrames
            concat_frames = []
            for key, path in paths.items():
                frame = read_frame(path)
                frame["root"] = key
                concat_frames.append(frame)

            combined_df = pd.concat(concat_frames, ignore_index=True)
            # Save the combined DataFrame to a file
            write_frames(combined_df, output_path, mode=mode)

        spec = CallSpec(
            path=("tables", "concat_tables"),
            kwargs={
                "paths": paths,
                "output_path": output_path,
                "mode": mode,
                "params": params,
            },
        ).render()
        return Runnable(__run, display=spec)

    @element
    def oldmerge(
        self,
        element: Element,
        annotation: Element,
        *,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        mode: str = "both",
        params: Params | None = None,
    ) -> Element:

        root = root or element.tag.root
        tag = from_prior(
            element.tag,
            root=root,
            state=State.ANNOTATED,
            method=Method.TABLES,
        )
        outdir = Path(outdir or element.file.parent)
        filename = filename or tag.default_output
        output_file = Path(outdir) / filename
        artifacts = Tables.artifacts_for_mode(output_file, mode=mode)

        run = self.old_merge_tables(
            element.file,
            annotation.file,
            left_name=element.tag.root,
            right_name=annotation.tag.root,
            output_path=output_file,
            mode=mode,
            params=params,
        )

        key, name = generate_element_key_name(tag, "Tables")

        return Element(
            key,
            run,
            tag=tag,
            determinants=(str(params),) if params else None,
            inputs=(element.file, annotation.file),
            artifacts=artifacts,
            pres=(element, annotation),
            name=name,
        )

    def old_merge_tables(
        self,
        left: Path,
        right: Path,
        *,
        left_name: str | None = None,
        right_name: str | None = None,
        output_path: Path,
        mode: str = "both",
        params: Params | None = None,
    ) -> Runnable:
        params = Params(
            how="left", suffixes=(f" ({left_name})", f" ({right_name})")
        ).update(params)

        def __run():
            # Concatenate the DataFrames
            left_df = read_frame(left)
            right_df = read_frame(right)
            combined_df = pd.merge(left_df, right_df, **params.to_dict())

            # Save the combined DataFrame to a file
            write_frames(combined_df, output_path, mode=mode)

        spec = CallSpec(
            path=("tables", "merge_tables"),
            kwargs={
                "left": left,
                "right": right,
                "output_path": output_path,
                "mode": mode,
                "params": params,
            },
        ).render()
        return Runnable(__run, display=spec)

    @element
    def merge(
        self,
        elements: Mapping[str, Element],
        *,
        filetags: Mapping[str, str] | None = None,
        root: str | None = None,
        how: Literal["left", "right", "outer", "inner"] = "left",
        on: str | Mapping[str, list[str]] | list[str] | None = None,
        select: list[str] | None = None,
        retains: bool = True,
        rename: Mapping[str, str] | None = None,
        order: list[str] | None = None,
        suffixmap: Mapping[str, str] | None = None,
        fillna: Any | None = 0,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        mode: str = "both",
        params: Params | None = None,
    ) -> Element:
        order = order or list(elements.keys())
        first_element = elements[order[0]] if order else next(iter(elements.values()))
        root = root or first_element.tag.root
        tag = from_prior(
            first_element.tag,
            tag,
            root=root,
            state=State.MERGED,
            method=Method.TABLES,
            param=how,
        )
        outdir = Path(outdir or first_element.file.parent)
        filename = filename or tag.default_output
        output_file = Path(outdir) / filename
        artifacts = Tables.artifacts_for_mode(output_file, mode=mode)
        if filetags:
            paths = {
                key: getattr(elements[key], filetags[key]) for key in elements.keys()
            }
        else:
            paths = {key: element.file for key, element in elements.items()}
        run = self.merge_tables(
            paths,
            output_path=output_file,
            how=how,
            on=on,
            select=select,
            retains=retains,
            order=order,
            rename=rename,
            suffixmap=suffixmap,
            fillna=fillna,
            mode=mode,
            params=params,
        )

        key, name = generate_element_key_name(tag, "Tables", subcommand="merge")

        return Element(
            key,
            run,
            tag=tag,
            determinants=(str(params),) if params else None,
            inputs=tuple(paths.values()),
            artifacts=artifacts,
            pres=tuple(elements.values()),
            name=name,
        )

    def merge_tables(
        self,
        paths: Mapping[str, Path],
        output_path: Path,
        *,
        how: Literal["left", "right", "outer", "inner"] = "left",
        on: str | Mapping[str, list[str]] | list[str] | None = None,
        select: list[str] | None = None,
        retains: bool = True,
        rename: Mapping[str, str] | None = None,
        order: list[str] | None = None,
        suffixmap: Mapping[str, str] | None = None,
        fillna: Any | None = 0,
        mode: str = "both",
        params: Params | None = None,
    ):

        @depends(drop_retain_rename)
        def __prep_frame(
            path: Path, merge_on: list[str] | None = None, suffix: str | None = None
        ) -> DataFrame:
            df = read_frame(path)
            selected = []
            if select:
                if retains:
                    selected = merge_on.copy() if merge_on else []
                    selected += select
                else:
                    if merge_on and set(merge_on).intersection(select):
                        raise ValueError(
                            f"Cannot drop columns {select} that are in 'on' columns {merge_on}"
                        )
                    selected = select

            morph = drop_retain_rename(
                select=selected,
                retains=retains,
                rename=rename,
            )
            df = morph(df)
            if suffix:
                # ensure the columns in first frames are renamed to include the key suffix
                rename_me = {
                    col: f"{col}{suffix}"
                    for col in df.columns
                    if col not in (merge_on or [])
                }
                print("renaming columns", rename_me)
                df = df.rename(columns=rename_me)
            print(df, merge_on)
            return df

        def __run(params=params):
            # Concatenate the DataFrames
            if isinstance(on, str):
                on_dict = {on: [on]}
            elif isinstance(on, dict):
                on_dict = on
            elif isinstance(on, list):
                on_dict = dict.fromkeys(paths.keys(), on)
            else:
                on_dict = {}
            first_key = order[0] if order else next(iter(paths.keys()))
            first_path = paths[first_key]
            # prepare dataframes to the relevant columns
            left_df = __prep_frame(
                first_path,
                merge_on=on_dict.get(first_key),
                suffix=suffixmap.get(first_key, None) if suffixmap else None,
            )

            print("+++++++++++++++")
            print("all the frames to join")
            for pp in paths.values():
                df = read_frame(pp)
                print(df.head(2))
            # iteratively merge the rest of the dataframes
            for key in order[1:] if order else list(paths.keys())[1:]:
                right_path = paths[key]
                right_df = __prep_frame(
                    right_path,
                    merge_on=on_dict.get(key),
                    suffix=suffixmap.get(key, None) if suffixmap else None,
                )
                if suffixmap:
                    suffixes = (
                        suffixmap.get(first_key, first_key),
                        suffixmap.get(key, key),
                    )
                else:
                    suffixes = (first_key, key)
                print(right_df.head(), "right_df", right_df.shape)
                left_on = on_dict.get(first_key)
                right_on = on_dict.get(key)
                validate_merge_columns(left_df, right_df, left_on, right_on)
                left_df = pd.merge(
                    left_df,
                    right_df,
                    how=how,
                    left_on=left_on,
                    right_on=right_on,
                    **params.to_dict() if params else {},
                )
                # drop the right_on column if it is the same as left_on
                if left_on != right_on and right_on is not None and left_on is not None:
                    left_df = left_df.drop(
                        columns=[col for col in right_on if col not in left_on]
                    )
            print("##############")
            # Save the combined DataFrame to a file
            columns_order = on_dict.get(first_key, []) + [
                col for col in left_df.columns if col not in on_dict.get(first_key, [])
            ]
            left_df = left_df[columns_order]
            left_df = left_df.fillna(fillna) if fillna is not None else left_df
            left_df = left_df.rename(columns=rename) if rename else left_df
            write_frames(left_df, output_path, mode=mode)

        spec = CallSpec(
            path=("tables", "merge_tables"),
            kwargs={
                "paths": paths,
                "output_path": output_path,
                "how": how,
                "on": on,
                "retained": retains,
                "order": order,
                "mode": mode,
                "params": params,
                "rename": rename,
            },
        ).render()
        return Runnable(__run, cmd=spec)

    @element
    def transform(
        self,
        element: TableElement,
        morph: Callable[[DataFrame], DataFrame],
        *,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        mode: str = "both",
        params: Params | None = None,
    ) -> Element:
        tag = from_prior(
            element.tag,
            tag,
            root=root,
            state=State.TRANSFORMED,
            method=Method.TABLES,
        )
        outdir = Path(outdir or element.file.parent)
        filename = filename or tag.default_output
        output_file = Path(outdir) / filename
        artifacts = Tables.artifacts_for_mode(output_file, mode=mode)

        def __run():
            df = read_frame(element.file)
            transformed_df = morph(df)
            write_frames(transformed_df, output_file, mode=mode)

        key, name = generate_element_key_name(tag, "Tables")

        return Element(
            key,
            __run,
            tag=tag,
            determinants=(str(params),) if params else None,
            inputs=(element.file,),
            artifacts=artifacts,
            pres=(element,),
            name=name,
        )

    @element
    def filter2annotate(
        self,
        element: Element,
        annotation: Element,
        *,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        mode: str = "both",
        params: Params | None = None,
        drop: list[str] | None = None,
        sortby: list[str] | None = None,
        ascending: list[bool] | bool = True,
        dtypes: dict[str, str] | None = None,
    ) -> Element:

        root = root or element.tag.root
        tag = from_prior(
            element.tag,
            tag,
            root=root,
            state=State.ANNOTATED,
            method=Method.TABLES,
        )
        outdir = Path(outdir or element.file.parent)
        filename = filename or tag.default_output
        output_file = Path(outdir) / filename
        artifacts = Tables.artifacts_for_mode(output_file, mode=mode)

        run = self.filter_to_annotate(
            annotation.file,
            element.file,
            left_name=annotation.tag.root,
            right_name=element.tag.root,
            output_path=output_file,
            mode=mode,
            params=params,
            drop=drop,
            sortby=sortby,
            ascending=ascending,
        )

        key, name = generate_element_key_name(tag, "Tables")

        return Element(
            key,
            run,
            tag=tag,
            determinants=(str(params),) if params else None,
            inputs=(element.file, annotation.file),
            artifacts=artifacts,
            pres=(element, annotation),
            name=name,
        )

    def filter_to_annotate(
        self,
        annotation: Path,
        right: Path,
        *,
        left_name: str | None = None,
        right_name: str | None = None,
        drop: list[str] | None = None,
        sortby: list[str] | None = None,
        ascending: list[bool] | bool = True,
        output_path: Path,
        mode: str = "both",
        params: Params | None = None,
    ) -> Runnable:
        params = Params(
            how="left", suffixes=(f" ({left_name})", f" ({right_name})")
        ).update(params)

        def __run():
            # Concatenate the DataFrames
            annotation_df = read_frame(annotation)
            right_df = read_frame(right)
            combined_df = pd.merge(annotation_df, right_df, **params.to_dict())
            dropme = drop or []
            if hasattr(params, "right_on"):
                col = (
                    params.left_on
                    if params.left_on in combined_df.columns
                    else params.left_on + params.suffixes[0]
                )
                if col in combined_df.columns:
                    dropme.append(col)
            if dropme:
                combined_df = combined_df.drop(columns=dropme)
            # Save the combined DataFrame to a file
            if sortby:
                combined_df = combined_df.sort_values(by=sortby, ascending=ascending)
            for col in combined_df.columns:
                if col.startswith("Count"):
                    combined_df[col] = combined_df[col].fillna(0)
                    combined_df[col] = combined_df[col].astype("int64")
            write_frames(combined_df, output_path, mode=mode)

        spec = CallSpec(
            path=("tables", "filter_to_annotate"),
            kwargs={
                "annotation": annotation,
                "right": right,
                "output_path": output_path,
                "mode": mode,
                "params": params,
            },
        ).render()
        return Runnable(__run, display=spec)

    # @element
    # def melt(
    #     self,
    #     element: TableElement,
    #     *,
    #     root: str | None = None,
    #     tag: PartialElementTag | ElementTag | None = None,
    #     outdir: Path | str | None = None,
    #     filename: Path | str | None = None,
    #     mode: str = "parquet",
    #     params: Params | None = None,
    # ) -> Element:

    #     root = root or element.tag.root
    #     tag = from_prior(
    #         element.tag,
    #         tag,
    #         root=root,
    #         state=State.TRANSFORMED,
    #         method=Method.TABLES,
    #         param="long",
    #     )
    #     outdir = Path(outdir or element.file.parent)
    #     filename = filename or tag.default_output
    #     output_file = Path(outdir) / filename
    #     artifacts = Tables.artifacts_for_mode(output_file, mode=mode)

    #     run = self.long_format(
    #         element.file,
    #         output_path=output_file,
    #         mode=mode,
    #         params=params,
    #     )

    #     key, name = generate_element_key_name(tag, "Tables")

    #     return Element(
    #         key,
    #         run,
    #         tag=tag,
    #         determinants=(str(params),) if params else None,
    #         inputs=(element.file,),
    #         artifacts=artifacts,
    #         pres=(element,),
    #         name=name,
    #     )


def drop_retain_rename(
    select: list[str] | None = None,
    retains: bool = True,
    rename: Mapping[str, str] | None = None,
) -> Callable[[pd.DataFrame], pd.DataFrame]:
    """
    Returns a function that selects or drops columns based on the provided list.
    If `retains` is True, only the columns in the `select` list are retained.
    If `retains` is False, the columns in the `select` list are dropped.
    If a rename mapping is provided, it will also rename the columns according to the mapping.
    """

    def _drop_frqs(df: pd.DataFrame) -> pd.DataFrame:
        if rename:
            df = df.rename(columns=rename)
        if select:
            print("this is select", select)
            for col in df.columns:
                print(
                    col,
                    any(sub in col for sub in select),
                    select[0] in col,
                    f"{select[0]} in {col}",
                    select[1] in col,
                    f"{select[1]} in {col}",
                )
            selected_columns = [
                col for col in df.columns.tolist() if any(sub in col for sub in select)
            ]
            print(
                selected_columns,
                "selected_columns",
                "retains",
                retains,
                "rename",
                rename,
            )
            if retains:
                assert all(
                    col in df.columns for col in selected_columns
                ), "Some columns in 'select' are not in the DataFrame"
                df = df[[col for col in df.columns if col in selected_columns]]
            else:
                assert all(
                    col in df.columns for col in selected_columns
                ), "Some columns in 'select' are not in the DataFrame"
                df = df[[col for col in df.columns if col not in selected_columns]]
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
            f"Consider renaming these columns in one of the DataFrames before merging."
        )
