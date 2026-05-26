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
from typing import Mapping

from mmalignments.models.elements import element
from mmalignments.models.tags import ElementTag, Method, Stage, State, PartialElementTag
from mmalignments.models.tags import from_prior

from .r_external import RScriptExternal
from .table_element import TableElement

logger = logging.getLogger(__name__)

_R_SCRIPT_DIR = Path(__file__).parent.parent.parent / "r"

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
            r_source_file=(r_script_dir or _R_SCRIPT_DIR) / "deseq2.R",
            version=version,
            source="https://bioconductor.org/packages/DESeq2",
        )

    # ------------------------------------------------------------------
    # Variance-stabilising transformation (VST)
    # ------------------------------------------------------------------

    @element
    def vst(
        self,
        counts: TableElement,
        view: str = "raw",
        *,
        sample_columns: Sequence[str] | None = None,
        blind: bool = True,
        fit_type: str = "parametric",
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
        propagate_roles: str | list[str] | None = "annotation",
    ) -> TableElement:
        """Variance-stabilising transformation (VST) via DESeq2 rpy2."""
        role = "vst"
        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.DESEQ2,
            state=State.NORMAL,
            ext="tsv",
        )

        counts_df = counts.view(view)
        if sample_columns:
            counts_df = counts_df[list(sample_columns)]
        cols = list(sample_columns) if sample_columns else list(counts_df.columns)

        out_dir = Path(outdir) if outdir else counts.tsv.parent
        output_tsv = out_dir / (filename or tag.default_output)
        output_parquet = output_tsv.with_suffix(".parquet")

        key, name = generate_element_key_name(tag, self.version_name, "vst")
        propagate_roles = [propagate_roles] if isinstance(propagate_roles, str) else propagate_roles
        if propagate_roles:
            propagation_df = counts.roles(*propagate_roles)
        column_roles: dict[str, str] = {c: role for c in cols}
        column_roles.update(counts.propagate("annotation"))
        runner = self.call_vst(
            counts_df=counts_df,
            output_tsv=output_tsv,
            output_parquet=output_parquet,
            sample_columns=cols,
            blind=blind,
            fit_type=fit_type,
            propagation_df=propagation_df,
            params=params,
            cfg=cfg,
        )
        return TableElement(
            key=key,
            run=runner,
            tag=tag,
            tsv=output_tsv,
            parquet=output_parquet,
            column_roles=column_roles,
            determinants=(blind, fit_type, self.version_name),
            inputs=(counts.tsv,),
            pres=(counts,),
            name=name,
        )

    @rsubroutine
    def call_vst(
        self,
        counts_df: DataFrame,
        output_tsv: Path,
        output_parquet: Path | None,
        sample_columns: Sequence[str] | None = None,
        blind: bool = True,
        fit_type: str = "parametric",
        propagation_df: DataFrame | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ):
        def post(result: DataFrame) -> DataFrame:
            print("VST result", result)
            result = result.copy()
            print("Propagation df", propagation_df)
            if isinstance(propagation_df, DataFrame):
                result = propagation_df.merge(result, left_index=True, right_index=True, how="right")

            result.to_csv(output_tsv, sep="\t", index=False)
            if output_parquet is not None:
                result.to_parquet(output_parquet, index=False)
            return result

        payload = {
            "counts_df": counts_df,
            "sample_columns": list(sample_columns) if sample_columns else None,
            "blind": blind,
            "fitType": fit_type,
        }
        return "counts_to_vst", payload, [], [output_tsv, output_parquet], None, None, post

    # ------------------------------------------------------------------
    # Regularised log (rlog)
    # ------------------------------------------------------------------

    @element
    def rlog(
        self,
        counts: TableElement,
        view: str = "raw_counts",
        *,
        sample_columns: Sequence[str] | None = None,
        blind: bool = True,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> TableElement:
        """Regularised log transformation (rlog) via DESeq2 rpy2."""
        role = "rlog"
        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.DESEQ2,
            state=State.NORMAL,
            ext="tsv",
        )

        counts_df = counts.view(view)
        if sample_columns:
            counts_df = counts_df[list(sample_columns)]
        cols = list(sample_columns) if sample_columns else list(counts_df.columns)

        out_dir = Path(outdir) if outdir else counts.tsv.parent
        output_tsv = out_dir / (filename or tag.default_output)
        output_parquet = output_tsv.with_suffix(".parquet")

        key, name = generate_element_key_name(tag, self.version_name, "rlog")

        column_roles: dict[str, str] = {c: role for c in cols}
        column_roles.update(counts.propagate())

        runner = self.call_rlog(
            counts_df=counts_df,
            output_tsv=output_tsv,
            output_parquet=output_parquet,
            sample_columns=cols,
            blind=blind,
            propagation_df=counts_df,
            params=params,
            cfg=cfg,
        )
        return TableElement(
            key=key,
            run=runner,
            tag=tag,
            tsv=output_tsv,
            parquet=output_parquet,
            column_roles=column_roles,
            determinants=(blind, self.version_name),
            inputs=(counts.parquet,),
            pres=(counts,),
            name=name,
        )

    @rsubroutine
    def call_rlog(
        self,
        counts_df: DataFrame,
        output_tsv: Path,
        output_parquet: Path | None,
        sample_columns: Sequence[str] | None = None,
        blind: bool = True,
        propagation_df: DataFrame | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ):
        def post(result: DataFrame) -> DataFrame:
            result = result.copy()
            if propagation_df is not None:
                for col in propagation_df.columns:
                    if col not in result.columns:
                        result[col] = propagation_df[col].values
            result.to_csv(output_tsv, sep="\t", index=True)
            if output_parquet is not None:
                result.to_parquet(output_parquet, index=True)
            return result

        payload = {
            "counts_df": counts_df,
            "sample_columns": list(sample_columns) if sample_columns else None,
            "blind": blind,
        }
        return "counts_to_rlog", payload, [], [output_tsv, output_parquet], None, None, post

    # ------------------------------------------------------------------
    # Unpaired two-group comparison
    # ------------------------------------------------------------------

    @element
    def unpaired(
        self,
        counts: TableElement,
        *,
        condition_a: str,
        condition_b: str,
        model_conditions: Mapping[str, Sequence[str]],
        comparison_name: str | None = None,
        column_map: Mapping[str, str] | None = None,
        view: str = "raw_counts",
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> TableElement:
        """Two-group Wald-test via DESeq2 rpy2."""
        comparison_name = comparison_name or f"{condition_a}_vs_{condition_b}"

        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.DIFF,
            method=Method.DESEQ2,
            state=State.DIFF,
            ext="tsv",
        )

        counts_df = counts.view(view)
        out_dir = Path(outdir) if outdir else counts.tsv.parent / "deseq2"
        stem = f"{counts.name}.deseq2.{comparison_name}"
        output_tsv = out_dir / f"{stem}.tsv"
        output_parquet = out_dir / f"{stem}.parquet"

        col_map = dict(column_map) if column_map else {
            k: v.format(comparison=comparison_name)
            for k, v in _DEFAULT_COLUMN_MAP_UNPAIRED.items()
        }

        key, name = generate_element_key_name(tag, self.version_name, comparison_name)

        all_sample_cols = [c for cols in model_conditions.values() for c in cols]
        column_roles: dict[str, str] = {c: "diff" for c in col_map.values()}
        column_roles.update(counts.propagate())

        runner = self.call_unpaired(
            counts_df=counts_df,
            condition_a=condition_a,
            condition_b=condition_b,
            model_conditions={k: list(v) for k, v in model_conditions.items()},
            column_map=col_map,
            output_tsv=output_tsv,
            output_parquet=output_parquet,
            propagation_df=counts_df,
            params=params,
            cfg=cfg,
        )
        return TableElement(
            key=key,
            run=runner,
            tag=tag,
            tsv=output_tsv,
            parquet=output_parquet,
            column_roles=column_roles,
            determinants=(comparison_name, condition_a, condition_b, self.version_name),
            inputs=(counts.parquet,),
            pres=(counts,),
            name=name,
        )

    @rsubroutine
    def call_unpaired(
        self,
        counts_df: DataFrame,
        condition_a: str,
        condition_b: str,
        model_conditions: dict[str, list[str]],
        column_map: dict[str, str],
        output_tsv: Path,
        output_parquet: Path | None,
        propagation_df: DataFrame | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ):
        def post(result: DataFrame) -> DataFrame:
            result = result.copy()
            if propagation_df is not None:
                for col in propagation_df.columns:
                    if col not in result.columns:
                        try:
                            result[col] = propagation_df.loc[result.index, col]
                        except Exception:
                            pass
            result.to_csv(output_tsv, sep="\t", index=True)
            if output_parquet is not None:
                result.to_parquet(output_parquet, index=True)
            return result

        payload = {
            "counts_df": counts_df,
            "condition_a": condition_a,
            "condition_b": condition_b,
            "model_conditions": model_conditions,
            "column_map": column_map,
        }
        return "deseq2_unpaired", payload, [], [output_tsv, output_parquet], None, None, post

    # ------------------------------------------------------------------
    # Time-series / multi-factor LRT
    # ------------------------------------------------------------------

    @element
    def timeseries(
        self,
        counts: TableElement,
        *,
        sample_columns: Sequence[str],
        factors: Mapping[str, Sequence[str]],
        formula: str,
        reduced: str,
        comparison_name: str | None = None,
        column_map: Mapping[str, str] | None = None,
        view: str = "raw_counts",
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> TableElement:
        """Multi-factor LRT via DESeq2 rpy2."""
        comparison_name = comparison_name or "timeseries"

        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.DIFF,
            method=Method.DESEQ2,
            state=State.DIFF,
            ext="tsv",
        )

        counts_df = counts.view(view)
        out_dir = Path(outdir) if outdir else counts.tsv.parent / "deseq2"
        stem = f"{counts.name}.deseq2.{comparison_name}"
        output_tsv = out_dir / f"{stem}.tsv"
        output_parquet = out_dir / f"{stem}.parquet"

        col_map = dict(column_map) if column_map else {
            k: v.format(comparison=comparison_name)
            for k, v in _DEFAULT_COLUMN_MAP_TIMESERIES.items()
        }

        key, name = generate_element_key_name(tag, self.version_name, comparison_name)

        column_roles: dict[str, str] = {c: "diff" for c in col_map.values()}
        column_roles.update(counts.propagate())

        runner = self.call_timeseries(
            counts_df=counts_df,
            sample_columns=list(sample_columns),
            factors={k: list(v) for k, v in factors.items()},
            formula=formula,
            reduced=reduced,
            column_map=col_map,
            output_tsv=output_tsv,
            output_parquet=output_parquet,
            propagation_df=counts_df,
            params=params,
            cfg=cfg,
        )
        return TableElement(
            key=key,
            run=runner,
            tag=tag,
            tsv=output_tsv,
            parquet=output_parquet,
            column_roles=column_roles,
            determinants=(comparison_name, formula, reduced, self.version_name),
            inputs=(counts.parquet,),
            pres=(counts,),
            name=name,
        )

    @rsubroutine
    def call_timeseries(
        self,
        counts_df: DataFrame,
        sample_columns: list[str],
        factors: dict[str, list[str]],
        formula: str,
        reduced: str,
        column_map: dict[str, str],
        output_tsv: Path,
        output_parquet: Path | None,
        propagation_df: DataFrame | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ):
        def post(result: DataFrame) -> DataFrame:
            result = result.copy()
            if propagation_df is not None:
                for col in propagation_df.columns:
                    if col not in result.columns:
                        try:
                            result[col] = propagation_df.loc[result.index, col]
                        except Exception:
                            pass
            result.to_csv(output_tsv, sep="\t", index=True)
            if output_parquet is not None:
                result.to_parquet(output_parquet, index=True)
            return result

        payload = {
            "counts_df": counts_df,
            "sample_columns": sample_columns,
            "factors": factors,
            "formula": formula,
            "reduced": reduced,
            "column_map": column_map,
        }
        return "deseq2_timeseries", payload, [], [output_tsv, output_parquet], None, None, post




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
    ...     condition_to_columns={"KO": ["KO_1", "KO_2"], "WT": ["WT_1", "WT_2"]},
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
            Maps condition name → list of sample column names in the count matrix.
        comparison_name : str | None
            Human-readable label used as column suffix, e.g. ``"KO_vs_WT"``.
            Defaults to ``"{condition_a}_vs_{condition_b}"``.
        include_other_columns_for_variance : bool
            If True, samples from other conditions are included in the DESeq2
            model (for improved dispersion estimates) but not in the contrast.
        tag : PartialElementTag | ElementTag | None
            Override for the output element tag.
        outdir : Path | str | None
            Output directory. Defaults to ``counts.tsv.parent / "differential"``.

        Returns
        -------
        TableElement
            Element with DESeq2 result columns.  Run it to materialise the files.
        """
        comparison_name = comparison_name or f"{condition_a}_vs_{condition_b}"
        suffix = _suffix(comparison_name, self.name)
        col_map = _colnames_deseq2(suffix)

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
        tsv_path     = outdir / f"{stem}.tsv"
        parquet_path = outdir / f"{stem}.parquet"
        params_path  = outdir / f"{stem}.params.json"

        # All condition columns that go into the model
        model_conditions: dict[str, list[str]] = {
            condition_a: list(condition_to_columns[condition_a]),
            condition_b: list(condition_to_columns[condition_b]),
        }
        if include_other_columns_for_variance:
            for cond, cols in condition_to_columns.items():
                if cond not in (condition_a, condition_b):
                    model_conditions[cond] = list(cols)

        payload = {
            "mode":               "unpaired",
            "counts_parquet":     counts.parquet,
            "condition_a":        condition_a,
            "condition_b":        condition_b,
            "model_conditions":   model_conditions,
            "comparison_name":    comparison_name,
            "column_map":         col_map,
            "output_tsv":         tsv_path,
            "output_parquet":     parquet_path,
        }
        self.write_params_json(params_path, payload)

        runner = self.run_r_script(
            script_name="deseq2.R",
            params_path=params_path,
            output_files=[tsv_path, parquet_path],
        )

        return TableElement(
            key=f"{counts.name}.deseq2_unpaired.{comparison_name}.{self.version_name}",
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

    @element
    def timeseries(
        self,
        counts: TableElement,
        *,
        sample_columns: Collection[str],
        factors: Mapping[str, list[str]],
        formula: str,
        reduced: str,
        comparison_name: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
    ) -> TableElement:
        """LRT-based DESeq2 analysis for time-series or multi-factor designs.

        Parameters
        ----------
        counts : TableElement
            Raw count matrix element.
        sample_columns : Collection[str]
            Ordered list of sample column names that appear in the count matrix.
        factors : Mapping[str, list[str]]
            Factor metadata passed to ``colData``.  Each key is a factor name
            (e.g. ``"condition"``, ``"time"``); values are per-sample labels in
            the same order as *sample_columns*.
        formula : str
            Full model design formula, e.g. ``"~ condition + time"``.
        reduced : str
            Reduced model formula for LRT, e.g. ``"~ condition"``.
        comparison_name : str | None
            Label used as column suffix and in file names.
        tag, outdir
            As in :meth:`unpaired`.

        Returns
        -------
        TableElement
        """
        comparison_name = comparison_name or "timeseries"
        suffix = _suffix(comparison_name, self.name)
        col_map = _colnames_deseq2(suffix)

        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.DESEQ2,
            state=State.DIFF,
            ext="tsv",
        )

        outdir = Path(outdir or counts.tsv.parent / "differential")
        stem = f"{counts.name}.deseq2_timeseries.{comparison_name}"
        tsv_path     = outdir / f"{stem}.tsv"
        parquet_path = outdir / f"{stem}.parquet"
        params_path  = outdir / f"{stem}.params.json"

        payload = {
            "mode":            "timeseries",
            "counts_parquet":  counts.parquet,
            "sample_columns":  list(sample_columns),
            "factors":         {k: list(v) for k, v in factors.items()},
            "formula":         formula,
            "reduced":         reduced,
            "comparison_name": comparison_name,
            "column_map":      col_map,
            "output_tsv":      tsv_path,
            "output_parquet":  parquet_path,
        }
        self.write_params_json(params_path, payload)

        runner = self.run_r_script(
            script_name="deseq2.R",
            params_path=params_path,
            output_files=[tsv_path, parquet_path],
        )

        return TableElement(
            key=f"{counts.name}.deseq2_timeseries.{comparison_name}.{self.version_name}",
            run=runner,
            tag=tag,
            tsv=tsv_path,
            parquet=parquet_path,
            determinants=(comparison_name, formula, reduced, self.version_name),
            inputs=(counts.parquet,),
            pres=(counts,),
            name=f"{counts.name}.deseq2_timeseries_{comparison_name}",
        )