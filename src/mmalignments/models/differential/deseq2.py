"""DESeq2 differential expression analysis — pipeline Element interface.

Each public method decorated with ``@element`` accepts a :class:`TableElement`
(count matrix) and returns a new :class:`TableElement` with DESeq2 results
columns appended (log2FC, p, FDR, …).

The actual computation is delegated to ``Rscript`` via
:meth:`RScriptExternal.run_r_script`.  Parameters are passed as a JSON file so
no rpy2 dependency is needed at the Python level.

Column naming convention
------------------------
All result columns carry a suffix derived from the comparison name so that
multiple comparisons can coexist in the same downstream table:

    ``log2FC (KO_vs_WT)``
    ``FDR (KO_vs_WT)``
    …
"""

from __future__ import annotations

import logging
from collections.abc import Collection
from pathlib import Path
from typing import Callable, Mapping, Sequence

from pandas import DataFrame

from mmalignments.core.annotations import _VIEWS, View
from mmalignments.models.artifacts import ArtifactSet, OutputSpec, TableArtifact
from mmalignments.models.elements import (
    Element,
    TableElement,
    element,
    generate_element_key_name,
)
from mmalignments.models.externals import (
    ExternalRunConfig,
)
from mmalignments.models.parameters import Params
from mmalignments.models.tables.frames import (
    Tables,
    TableSource,
    as_table_source,
    divide,
)
from mmalignments.models.tags import (
    ElementTag,
    Method,
    PartialElementTag,
    Stage,
    State,
    from_prior,
)
from mmalignments.r import (
    R_SCRIPT_DIR,
    RScriptExternal,
    RScriptInternal,
    RSubroutineIn,
    rsubroutine,
)
from mmalignments.services.io import read_frame, write_frame, write_frames

logger = logging.getLogger(__name__)


# Default column renaming: DESeq2 R column -> pipeline column template.
# The comparison name is inserted at format time.
_DEFAULT_COLUMN_MAP_UNPAIRED: dict[str, str] = {
    "log2FoldChange": "log2FC ({comparison})",
    "pvalue": "p ({comparison})",
    "padj": "FDR ({comparison})",
    "baseMean": "baseMean ({comparison})",
    "lfcSE": "lfcSE ({comparison})",
    "stat": "stat ({comparison})",
}

_DEFAULT_COLUMN_MAP_TIMESERIES: dict[str, str] = {
    "baseMean": "baseMean ({comparison})",
    "stat": "stat ({comparison})",
    "pvalue": "p ({comparison})",
    "padj": "FDR ({comparison})",
}


class DESeq2In(RScriptInternal):
    """Pipeline interface for DESeq2 using rpy2 (in-process, no subprocess).

    R functions are called directly from ``deseq2.R`` via rpy2.
    File I/O (TSV / Parquet) and annotation propagation are handled in Python
    post-callbacks.

    Examples
    --------
    >>> deseq2 = DESeq2In()
    >>> vst_el = deseq2.vst(counts_element, view="raw_counts")
    >>> rlog_el = deseq2.rlog(counts_element, view="raw_counts")
    >>> result  = deseq2.unpaired(
    ...     counts_element,
    ...     condition_a="KO",
    ...     condition_b="WT",
    ...     model_conditions={"KO": ["KO_1", "KO_2"], "WT": ["WT_1", "WT_2"]},
    ... )
    """

    def __init__(
        self,
        r_script_dir: Path | None = None,
        version: str | None = None,
    ) -> None:
        super().__init__(
            name="DESeq2",
            r_source_file=(r_script_dir or R_SCRIPT_DIR) / "deseq2.R",
            version=version,
            source="https://bioconductor.org/packages/DESeq2",
        )
        self.default_output_spec = OutputSpec(
            ext="parquet", additional_extensions=["tsv"]
        )

    @property
    def default_dir(self) -> Path:
        return Path("results/deseq2")

    ############################################################################
    # Helper
    ############################################################################

    def merge_rename_write(
        self,
        outputs: list[Path],
        other: DataFrame | None = None,
        rename: dict[str, str] | None = None,
    ) -> Callable[[DataFrame], DataFrame]:
        """
        Return a post-processing function that merges the result with other
        and writes to disk.
        """

        def post(result: DataFrame) -> DataFrame:
            if rename:
                result = result.rename(columns=rename)
            return merge_and_write(
                result=result,
                outputs=outputs,
                other_df=other,
            )

        return post

    # ------------------------------------------------------------------
    # Unpaired two-group comparison
    # ------------------------------------------------------------------
    @element
    def unpaired(
        self,
        counts: TableElement,
        condition_a: str,
        condition_b: str,
        model_conditions: Mapping[str, list[str]],
        *,
        comparison_name: str | None = None,
        # column_schema: ColumnSchema | None = None,
        view: View | None = None,
        use_all_samples_for_variance: bool = True,
        propagation: bool = False,  # if the result is appended or not, do this differently # noqa: E501
        tag: PartialElementTag | ElementTag | None = None,
        output_spec: OutputSpec | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Two-group Wald-test via DESeq2 rpy2."""
        comparison_name = comparison_name or f"{condition_a}_vs_{condition_b}"
        output_spec = output_spec or self.default_output_spec
        view = view or _VIEWS["escd"]
        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.DIFF,
            method=Method.DESEQ2,
            state=State.DIFF,
        )
        infile = counts.primary.resolve()

        col_map = {
            k: v.format(comparison=comparison_name)
            for k, v in _DEFAULT_COLUMN_MAP_UNPAIRED.items()
        }

        # column_schema = column_schema or ColumnSchema.build(col_map.values())
        artifacts, output = ArtifactSet.generate_file_artifacts(
            tag=tag,
            infile=infile,
            spec=output_spec,
            # column_schema=column_schema,
            default_dir=self.default_dir,
        )
        artifacts["size_factors"] = TableArtifact(
            path=output.with_name("size_factors"),
        )
        additional_filter = None
        if not use_all_samples_for_variance:
            additional_filter = filter_to_condition_sample(
                model_conditions, condition_a, condition_b
            )
        key, name = generate_element_key_name(
            tag, self.version_name, comparison_name
        )  # noqa: E501
        runner = self.call_unpaired(
            source=as_table_source(counts),
            condition_a=condition_a,
            condition_b=condition_b,
            model_conditions={k: list(v) for k, v in model_conditions.items()},
            column_map=col_map,
            output=output,
            view=view,
            propagation=propagation,
            additional_filter=additional_filter,
            cfg=cfg,
        )
        return Element(
            key=key,
            run=runner,
            tag=tag,
            artifacts=artifacts,
            determinants=(
                comparison_name,
                condition_a,
                condition_b,
                self.version_name,
                str(view),
            ),
            inputs=(infile,),
            pres=(counts,),
            name=name,
        )

    @rsubroutine
    def call_unpaired(
        self,
        source: TableSource,
        condition_a: str,
        condition_b: str,
        model_conditions: dict[str, list[str]],
        column_map: dict[str, str],
        output: Path | list[Path],
        view: View | None = None,
        propagation: bool = False,
        additional_filter: Callable[[DataFrame], DataFrame] | None = None,
        cfg: ExternalRunConfig | None = None,
        # params: Params | None = None,
    ) -> RSubroutineIn:
        outputs = output if isinstance(output, list) else [output]

        def make_context():
            df = source.view(view)
            if additional_filter:
                df = additional_filter(df)
            counts_df = df[
                model_conditions[condition_a] + model_conditions[condition_b]
            ]
            other_df = (
                df[[c for c in df.columns if c not in counts_df.columns]]
                if propagation
                else None
            )

            return {
                "payload": {
                    "counts_df": counts_df,
                    "condition_a": condition_a,
                    "condition_b": condition_b,
                    "model_conditions": model_conditions,
                    "column_map": column_map,
                },
                "context": {"other_df": other_df, "outputs": outputs},
            }

        def post(result, context):
            results = result["results"]
            size_factors = result["size_factors"]
            other_df = context["other_df"]
            outputs = context["outputs"]
            fn = self.merge_rename_write(outputs=outputs, other=other_df)
            write_frame(
                size_factors, [outputs[0].with_name("size_factors")]
            )  # noqa: E501
            fn(results)

        cfg = cfg or ExternalRunConfig(threads=1)
        return (
            "deseq2_unpaired",
            make_context,
            [source.resolve()],
            outputs,
            None,
            None,
            post,
        )

    # ------------------------------------------------------------------
    # Variance-stabilising transformation (VST)
    # ------------------------------------------------------------------

    @element
    def vst(
        self,
        counts: Element,
        *,
        blind: bool = True,
        fit_type: str = "parametric",
        sample_columns: Sequence[str] | None = None,
        view: View | None = None,
        propagation: bool = False,
        tag: PartialElementTag | ElementTag | None = None,
        output_spec: OutputSpec | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Variance-stabilising transformation (VST) via DESeq2 rpy2."""
        if sample_columns is None and view is None:
            raise ValueError(
                "Either sample_columns or view must be provided to select the counts columns."  # noqa: E501
            )
        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.DESEQ2,
            state=State.NORMAL,
            ext="parquet",
            param="vst",
        )
        infile = counts.primary.resolve()
        # columns_used = counts.primary.column_schema.select(view)
        # new_schema, rename = counts.primary.column_schema.derive(
        #     columns_used,
        #     add_metric="VST",
        # )
        artifacts, output = ArtifactSet.generate_file_artifacts(
            tag=tag,
            infile=infile,
            spec=output_spec or self.default_output_spec,
            # column_schema=new_schema,
        )

        # df_call = frame_callable(counts.primary, view=view)

        # if sample_columns:
        #     counts_df = counts_df[list(sample_columns)]
        # cols = list(sample_columns) if sample_columns else list(counts_df.columns)    # noqa: E501

        key, name = generate_element_key_name(tag, self.version_name, "vst")
        runner = self.call_vst(
            source=as_table_source(counts),
            output=output,
            blind=blind,
            fit_type=fit_type,
            sample_columns=sample_columns,
            view=view,
            # rename=rename,
            propagation=counts.primary.load if propagation else None,
            # params=params,
            cfg=cfg,
        )
        return Element(
            key=key,
            run=runner,
            tag=tag,
            determinants=(str(blind), fit_type, self.version_name, str(view)),
            inputs=(infile,),
            artifacts=artifacts,
            pres=(counts,),
            name=name,
        )

    @rsubroutine
    def call_vst(
        self,
        source: TableSource,
        output: Path | list[Path],
        *,
        blind: bool = True,
        fit_type: str = "parametric",
        sample_columns: Sequence[str] | None = None,
        view: View | None = None,
        rename: dict[str, str] | None = None,
        propagation: bool = False,  # if the result is appended or not
        # params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ):

        outputs = output if isinstance(output, list) else [output]
        cfg = cfg or ExternalRunConfig(threads=1)

        def make_context():
            count_df = (
                source.view(view)[sample_columns]
                if sample_columns
                else source.view(view)
            )
            other_df = (
                source.frame[
                    [
                        c
                        for c in source.frame.columns
                        if c not in count_df.columns  # noqa: E501
                    ]  # noqa: E501
                ]
                if propagation
                else None
            )
            return {
                "payload": {
                    "counts_df": count_df,
                    "sample_columns": (
                        list(sample_columns) if sample_columns else None
                    ),  # noqa: E501
                    "blind": blind,
                    "fitType": fit_type,
                },
                "context": {"other_df": other_df, "rename": rename},
            }

        def post(result: DataFrame, context: dict) -> DataFrame:
            rename = context["context"].get("rename", None)
            other_df = context["context"].get("other_df", None)
            if rename:
                result = result.rename(columns=rename)
            return merge_and_write(
                result=result,
                outputs=outputs,
                other_df=other_df,
            )

        return (
            "counts_to_vst",
            make_context,
            [],
            outputs,
            None,
            None,
            post,
        )

    def normalize(
        source: Element,
        deseq: Element,
        *,
        root: str | None = None,
        index_column: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        views: Mapping[str, View] | None = None,
        output_spec: OutputSpec | None = None,
        params: Params | None = None,
    ) -> Element:

        morph = divide(
            numerator=source.root, denominator=deseq.root, views=views
        )  # noqa: E501
        return Tables().combine(
            source,
            deseq,
            morph,
            root=root,
            index_column=index_column,
            tag=tag,
            views=views,
            output_spec=output_spec,
            params=params,
        )

    # ------------------------------------------------------------------
    # Regularised log (rlog)
    # ------------------------------------------------------------------

    # @element
    # def rlog(
    #     self,
    #     counts: TableElement,
    #     view: str = "raw_counts",
    #     *,
    #     sample_columns: Sequence[str] | None = None,
    #     blind: bool = True,
    #     tag: PartialElementTag | ElementTag | None = None,
    #     outdir: Path | str | None = None,
    #     filename: str | None = None,
    #     params: Params | None = None,
    #     cfg: ExternalRunConfig | None = None,
    # ) -> TableElement:
    #     """Regularised log transformation (rlog) via DESeq2 rpy2."""
    #     role = "rlog"
    #     tag = from_prior(
    #         counts.tag,
    #         tag,
    #         stage=Stage.ANALYSIS,
    #         method=Method.DESEQ2,
    #         state=State.NORMAL,
    #         ext="tsv",
    #     )

    #     counts_df = counts.view(view)
    #     if sample_columns:
    #         counts_df = counts_df[list(sample_columns)]
    #     cols = list(sample_columns) if sample_columns else list(counts_df.columns)    # noqa: E501

    #     out_dir = Path(outdir) if outdir else counts.tsv.parent
    #     output_tsv = out_dir / (filename or tag.default_output)
    #     output_parquet = output_tsv.with_suffix(".parquet")

    #     key, name = generate_element_key_name(tag, self.version_name, "rlog")

    #     column_Role: dict[str, str] = dict.fromkeys(cols, role)
    #     column_Role.update(counts.propagate())

    #     runner = self.call_rlog(
    #         counts_df=counts_df,
    #         output_tsv=output_tsv,
    #         output_parquet=output_parquet,
    #         sample_columns=cols,
    #         blind=blind,
    #         propagation_df=counts_df,
    #         params=params,
    #         cfg=cfg,
    #     )
    #     return TableElement(
    #         key=key,
    #         run=runner,
    #         tag=tag,
    #         tsv=output_tsv,
    #         parquet=output_parquet,
    #         column_Role=column_Role,
    #         determinants=(blind, self.version_name),
    #         inputs=(counts.parquet,),
    #         pres=(counts,),
    #         name=name,
    #     )

    # @rsubroutine
    # def call_rlog(
    #     self,
    #     counts_df: DataFrame,
    #     output_tsv: Path,
    #     output_parquet: Path | None,
    #     sample_columns: Sequence[str] | None = None,
    #     blind: bool = True,
    #     propagation_df: DataFrame | None = None,
    #     params: Params | None = None,
    #     cfg: ExternalRunConfig | None = None,
    # ):
    #     def post(result: DataFrame) -> DataFrame:
    #         return merge_and_write(
    #             result=result,
    #             full_df=full_df,
    #             outputs=outputs,
    #             propagation=propagation,
    #         )

    #     payload = {
    #         "counts_df": counts_df,
    #        "sample_columns": list(sample_columns) if sample_columns else None,
    #         "blind": blind,
    #     }
    #     return (
    #         "counts_to_rlog",
    #         payload,
    #         [],
    #         [output_tsv, output_parquet],
    #         None,
    #         None,
    #         post,
    #     )

    # ------------------------------------------------------------------
    # Time-series / multi-factor LRT
    # ------------------------------------------------------------------

    # @element
    # def timeseries(
    #     self,
    #     counts: TableElement,
    #     *,
    #     sample_columns: Sequence[str],
    #     factors: Mapping[str, Sequence[str]],
    #     formula: str,
    #     reduced: str,
    #     comparison_name: str | None = None,
    #     column_map: Mapping[str, str] | None = None,
    #     view: str = "raw_counts",
    #     tag: PartialElementTag | ElementTag | None = None,
    #     outdir: Path | str | None = None,
    #     params: Params | None = None,
    #     cfg: ExternalRunConfig | None = None,
    # ) -> TableElement:
    #     """Multi-factor LRT via DESeq2 rpy2."""
    #     comparison_name = comparison_name or "timeseries"

    #     tag = from_prior(
    #         counts.tag,
    #         tag,
    #         stage=Stage.DIFF,
    #         method=Method.DESEQ2,
    #         state=State.DIFF,
    #         ext="tsv",
    #     )

    #     counts_df = counts.view(view)
    #     out_dir = Path(outdir) if outdir else counts.tsv.parent / "deseq2"
    #     stem = f"{counts.name}.deseq2.{comparison_name}"
    #     output_tsv = out_dir / f"{stem}.tsv"
    #     output_parquet = out_dir / f"{stem}.parquet"

    #     col_map = (
    #         dict(column_map)
    #         if column_map
    #         else {
    #             k: v.format(comparison=comparison_name)
    #             for k, v in _DEFAULT_COLUMN_MAP_TIMESERIES.items()
    #         }
    #     )

    #     key, name = generate_element_key_name(tag, self.version_name, comparison_name)    # noqa: E501

    #     column_Role: dict[str, str] = dict.fromkeys(col_map.values(), "diff")
    #     column_Role.update(counts.propagate())

    #     runner = self.call_timeseries(
    #         counts_df=counts_df,
    #         sample_columns=list(sample_columns),
    #         factors={k: list(v) for k, v in factors.items()},
    #         formula=formula,
    #         reduced=reduced,
    #         column_map=col_map,
    #         output_tsv=output_tsv,
    #         output_parquet=output_parquet,
    #         propagation_df=counts_df,
    #         params=params,
    #         cfg=cfg,
    #     )
    #     return TableElement(
    #         key=key,
    #         run=runner,
    #         tag=tag,
    #         tsv=output_tsv,
    #         parquet=output_parquet,
    #         column_Role=column_Role,
    #         determinants=(comparison_name, formula, reduced, self.version_name),  # noqa: E501
    #         inputs=(counts.parquet,),
    #         pres=(counts,),
    #         name=name,
    #     )

    # @rsubroutine
    # def call_timeseries(
    #     self,
    #     counts_df: DataFrame,
    #     sample_columns: list[str],
    #     factors: dict[str, list[str]],
    #     formula: str,
    #     reduced: str,
    #     column_map: dict[str, str],
    #     output_tsv: Path,
    #     output_parquet: Path | None,
    #     propagation_df: DataFrame | None = None,
    #     params: Params | None = None,
    #     cfg: ExternalRunConfig | None = None,
    # ):
    #     def post(result: DataFrame) -> DataFrame:
    #         return merge_and_write(
    #             result=result,
    #             full_df=full_df,
    #             outputs=outputs,
    #             propagation=propagation,
    #         )

    #     payload = {
    #         "counts_df": counts_df,
    #         "sample_columns": sample_columns,
    #         "factors": factors,
    #         "formula": formula,
    #         "reduced": reduced,
    #         "column_map": column_map,
    #     }
    #     return (
    #         "deseq2_timeseries",
    #         payload,
    #         [],
    #         [output_tsv, output_parquet],
    #         None,
    #         None,
    #         post,
    #     )


################################################################################
# helper
################################################################################


def merge_and_write(
    result: DataFrame,
    outputs: list[Path],
    other_df: DataFrame | None = None,
) -> DataFrame:
    if other_df is not None:
        result = other_df.loc[result.index, result.columns] = result
    write_frames(result, outputs)
    return result


def counts_call(
    source: Path | TableArtifact,
    view: View | None,
    additional_filter: Callable[[DataFrame], DataFrame] | None = None,
) -> Callable:
    """
    Prepare a counts DataFrame for DESeq2 analysis.

    This function returns a callable that, when invoked, reads the counts data
    from the specified source, applies the additional filter, and returns both
    the filtered counts DataFrame and the full DataFrame for propagation.

    Parameters
    ----------
    source : Path | TableArtifact
        The source of the counts data. It can be a Path to a file or a
        TableArtifact.
    view : View | None
        The view to filter for in the TableArtifact. If None, no filtering is
        applied.
    additional_filter : Callable[[DataFrame], DataFrame]
        A function that takes a DataFrame and returns a filtered DataFrame. This
        is applied to the counts DataFrame after reading it from the source.

    Returns
    -------
    Callable
        A callable that returns a tuple of
        (filtered counts DataFrame, full DataFrame).
    """

    def __counts_call() -> tuple[DataFrame, DataFrame]:
        """
        Return the count DataFrame containing only relevant data and the full
        frame for propagation.
        """
        if isinstance(source, Path):
            counts_df = read_frame(source)
            full = counts_df.copy()
        elif isinstance(source, TableArtifact):
            counts_df = source.view(view)
            full = source.view(view)
        else:
            raise ValueError(
                f"Unsupported counts type: {type(source)}. "
                "Expected Path or TableArtifact."
            )
        if additional_filter is not None:
            counts_df = additional_filter(counts_df)
        return counts_df, full

    return __counts_call


def filter_to_condition_sample(
    model_conditions: Mapping[str, list[str]],
    condition_a: str,
    condition_b: str,  # noqa: E501
) -> Callable[[DataFrame], DataFrame]:
    """
    Filter the counts DataFrame to only include columns corresponding to the
    specified conditions.

    Parameters
    ----------
    model_conditions : Mapping[str, list[str]]
        Dictionary mapping conditions to the corresponding sample columns.
    condition_a : str
        Name of the first condition.
    condition_b : str
        Name of the second condition.

    Returns
    -------
    Callable[[DataFrame], DataFrame]
        Function that filters a counts DataFrame to only include columns for the
        specified conditions.
    """

    def __filter(counts_df: DataFrame) -> DataFrame:
        """
        Filter the counts DataFrame to only include columns corresponding to the
        specified conditions.
        """
        drop_cols = [
            c
            for c in counts_df.columns
            if c
            not in model_conditions[condition_a]
            + model_conditions[condition_b]  # noqa: E501
        ]
        counts_df = counts_df.drop(columns=drop_cols)
        return counts_df

    return __filter


class DESeq2(RScriptExternal):
    """Pipeline interface for DESeq2 differential expression analysis.

    Two analysis modes are available:

    * :meth:`unpaired`  — standard two-group contrast (``DESeq2Unpaired``
      from the old code).
    * :meth:`timeseries` — LRT-based multi-factor/time-series design.

    Examples
    --------
    >>> deseq = DESeq2()
    >>> result_elem = deseq.unpaired(
    ...     counts_element,
    ...     condition_a="KO",
    ...     condition_b="WT",
    ...     condition_to_columns={"KO": ["KO_1", "KO_2"], "WT": ["WT_1",...]},
    ...     comparison_name="KO_vs_WT",
    ...     outdir=Path("results/differential"),
    ... )
    """

    def __init__(
        self,
        r_script_dir: Path | None = None,
        version: str | None = None,
    ) -> None:
        super().__init__(
            name="DESeq2",
            r_script_dir=r_script_dir,
            version=version,
            source="https://bioconductor.org/packages/DESeq2",
        )

    # ------------------------------------------------------------------
    # Unpaired two-group contrast
    # ------------------------------------------------------------------

    @element
    def unpaired(
        self,
        counts: TableElement,
        *,
        condition_a: str,
        condition_b: str,
        condition_to_columns: Mapping[str, Collection[str]],
        comparison_name: str | None = None,
        include_other_columns_for_variance: bool = False,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
    ) -> TableElement:
        """Two-group unpaired DESeq2 contrast.

        Parameters
        ----------
        counts : TableElement
            Element whose ``parquet`` artefact is the raw count matrix
            (genes × samples, integer counts).
        condition_a : str
            Name of the treatment / test condition.
        condition_b : str
            Name of the reference condition.
        condition_to_columns : Mapping[str, Collection[str]]
            Maps condition name → list of sample column names in the count
            matrix.
        comparison_name : str | None
            Human-readable label used as column suffix, e.g. ``"KO_vs_WT"``.
            Defaults to ``"{condition_a}_vs_{condition_b}"``.
        include_other_columns_for_variance : bool
            If True, samples from other conditions are included in the DESeq2
            model (for improved dispersion estimates) but not in the contrast.
        tag : PartialElementTag | ElementTag | None
            Override for the output element tag.
        outdir : Path | str | None
            Output directory. Defaults to ``counts.tsv.parent/"differential"``.

        Returns
        -------
        TableElement
            Element with DESeq2 result columns. Run it to materialise the files.
        """
        comparison_name = comparison_name or f"{condition_a}_vs_{condition_b}"
        # suffix = _suffix(comparison_name, self.name)
        # col_map = _colnames_deseq2(suffix)

        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.DESEQ2,
            state=State.DIFF,
            ext="tsv",
        )

        outdir = Path(outdir or counts.tsv.parent / "differential")
        stem = f"{counts.name}.deseq2_unpaired.{comparison_name}"
        tsv_path = outdir / f"{stem}.tsv"
        parquet_path = outdir / f"{stem}.parquet"
        params_path = outdir / f"{stem}.params.json"

        # All condition columns that go into the model
        model_conditions: dict[str, list[str]] = {
            condition_a: list(condition_to_columns[condition_a]),
            condition_b: list(condition_to_columns[condition_b]),
        }
        if include_other_columns_for_variance:
            for cond, cols in condition_to_columns.items():
                if cond not in (condition_a, condition_b):
                    model_conditions[cond] = list(cols)
        col_map = {
            k: v.format(comparison=comparison_name)
            for k, v in _DEFAULT_COLUMN_MAP_UNPAIRED.items()
        }

        payload = {
            "mode": "unpaired",
            "counts_parquet": counts.parquet,
            "condition_a": condition_a,
            "condition_b": condition_b,
            "model_conditions": model_conditions,
            "comparison_name": comparison_name,
            "column_map": col_map,
            "output_tsv": tsv_path,
            "output_parquet": parquet_path,
        }
        self.write_params_json(params_path, payload)

        runner = self.run_r_script(
            script_name="deseq2.R",
            params_path=params_path,
            output_files=[tsv_path, parquet_path],
        )

        return TableElement(
            key=f"{counts.name}.deseq2_unpaired.{comparison_name}.{self.version_name}",  # noqa: E501
            run=runner,
            tag=tag,
            tsv=tsv_path,
            parquet=parquet_path,
            determinants=(
                comparison_name,
                condition_a,
                condition_b,
                include_other_columns_for_variance,
                self.version_name,
            ),
            inputs=(counts.parquet,),
            pres=(counts,),
            name=f"{counts.name}.deseq2_{comparison_name}",
        )

    # ------------------------------------------------------------------
    # Timeseries / LRT design
    # ------------------------------------------------------------------

    # @element
    # def timeseries(
    #     self,
    #     counts: TableElement,
    #     *,
    #     sample_columns: Collection[str],
    #     factors: Mapping[str, list[str]],
    #     formula: str,
    #     reduced: str,
    #     comparison_name: str | None = None,
    #     tag: PartialElementTag | ElementTag | None = None,
    #     outdir: Path | str | None = None,
    # ) -> TableElement:
    #     """LRT-based DESeq2 analysis for time-series or multi-factor designs.

    #     Parameters
    #     ----------
    #     counts : TableElement
    #         Raw count matrix element.
    #     sample_columns : Collection[str]
    #       Ordered list of sample column names that appear in the count matrix.
    #     factors : Mapping[str, list[str]]
    #         Factor metadata passed to ``colData``.  Each key is a factor name
    #       (e.g. ``"condition"``, ``"time"``); values are per-sample labels in
    #         the same order as *sample_columns*.
    #     formula : str
    #         Full model design formula, e.g. ``"~ condition + time"``.
    #     reduced : str
    #         Reduced model formula for LRT, e.g. ``"~ condition"``.
    #     comparison_name : str | None
    #         Label used as column suffix and in file names.
    #     tag, outdir
    #         As in :meth:`unpaired`.

    #     Returns
    #     -------
    #     TableElement
    #     """
    #     comparison_name = comparison_name or "timeseries"
    #     suffix = _suffix(comparison_name, self.name)
    #     col_map = _colnames_deseq2(suffix)

    #     tag = from_prior(
    #         counts.tag,
    #         tag,
    #         stage=Stage.ANALYSIS,
    #         method=Method.DESEQ2,
    #         state=State.DIFF,
    #         ext="tsv",
    #     )

    #     outdir = Path(outdir or counts.tsv.parent / "differential")
    #     stem = f"{counts.name}.deseq2_timeseries.{comparison_name}"
    #     tsv_path = outdir / f"{stem}.tsv"
    #     parquet_path = outdir / f"{stem}.parquet"
    #     params_path = outdir / f"{stem}.params.json"

    #     payload = {
    #         "mode": "timeseries",
    #         "counts_parquet": counts.parquet,
    #         "sample_columns": list(sample_columns),
    #         "factors": {k: list(v) for k, v in factors.items()},
    #         "formula": formula,
    #         "reduced": reduced,
    #         "comparison_name": comparison_name,
    #         "column_map": col_map,
    #         "output_tsv": tsv_path,
    #         "output_parquet": parquet_path,
    #     }
    #     self.write_params_json(params_path, payload)

    #     runner = self.run_r_script(
    #         script_name="deseq2.R",
    #         params_path=params_path,
    #         output_files=[tsv_path, parquet_path],
    #     )

    #     return TableElement(


# key=f"{counts.name}.deseq2_timeseries.{comparison_name}.{self.version_name}",
#         run=runner,
#         tag=tag,
#         tsv=tsv_path,
#         parquet=parquet_path,
#       determinants=(comparison_name, formula, reduced, self.version_name),
#         inputs=(counts.parquet,),
#         pres=(counts,),
#         name=f"{counts.name}.deseq2_timeseries_{comparison_name}",
#     )
