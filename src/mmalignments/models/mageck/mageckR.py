"""MAGeCKFlute scatter-view analysis — pipeline Element interface.

Each public method decorated with ``@element`` accepts a :class:`TableElement`
(MAGeCK MLE or RRA result table) and returns an :class:`Element` whose
artefacts are plot files (JPEG) and a selection TSV of hit genes.

The actual computation is delegated to ``Rscript mageck.R`` via
:meth:`RScriptExternal.run_r_script`.  All parameters are passed as a JSON
file, so no rpy2 dependency is needed at the Python level.

Available analyses
------------------
* :meth:`scatter_view` — ScatterView hit calling from two beta / lfc columns.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

from mmalignments.models.elements import Element, element
from mmalignments.models.tags import (
    Method,
    Stage,
    State,
)
from mmalignments.models.overlay import (
    ElementTag,
    PartialElementTag,
    from_prior,
)

from mmalignments.r.r_integration import RScriptExternal

logger = logging.getLogger(__name__)

_R_SCRIPT_DIR = Path(__file__).parent.parent.parent / "r"


class MAGeCKR(RScriptExternal):
    """Pipeline interface for MAGeCKFlute visualisation and hit calling.

    Examples
    --------
    >>> mageck = MAGeCKR()
    >>> scatter = mageck.scatter_view(
    ...     mle_result_element,
    ...     x_col="WT|beta",
    ...     y_col="KO|beta",
    ...     data_type="mle",
    ...     comparison_name="KO_vs_WT",
    ... )
    """

    def __init__(
        self,
        r_script_dir: Path | None = None,
        version: str | None = None,
    ) -> None:
        super().__init__(
            name="MAGeCK",
            r_script_dir=r_script_dir or _R_SCRIPT_DIR,
            version=version,
            source="https://bioconductor.org/packages/MAGeCKFlute",
        )

    # ------------------------------------------------------------------
    # ScatterView
    # ------------------------------------------------------------------

    @element
    def scatter_view(
        self,
        table: Element,
        *,
        x_col: str,
        y_col: str,
        data_type: str = "mle",
        gene_col: str = "id",
        comparison_name: str | None = None,
        normalize: bool = True,
        top: int = 10,
        groups: Sequence[str] | None = None,
        select: str | None = None,
        neg_effect_cutoff: float = -0.4,
        pos_effect_cutoff: float = 0.4,
        delta_cutoff_k: float = 2.0,
        filter_fdr_x: bool = False,
        filter_fdr_y: bool = False,
        filter_groups: bool = True,
        fdr_cutoff: float = 0.05,
        toplabels: bool = True,
        label_selected_only: bool = False,
        xlab: str | None = None,
        ylab: str | None = None,
        jpeg_width: int = 20,
        jpeg_height: int = 15,
        auto_cut_diag: float = 2.0,
        auto_cut_x: float = 2.0,
        auto_cut_y: float = 2.0,
        x_cut: tuple[float, float] | None = None,
        y_cut: tuple[float, float] | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
    ) -> Element:
        """Scatter-view hit calling via MAGeCKFlute.

        Reads a MAGeCK MLE or RRA result TSV, produces a ``ScatterView``
        JPEG, a QC linear-fit JPEG, a data TSV with enrichment annotations,
        and a TSV of selected hit genes.

        Parameters
        ----------
        table : Element
            Input element.  Must expose a ``tsv`` artefact pointing to the
            MAGeCK result table.
        x_col : str
            Column name for the X-axis (e.g. ``"WT|beta"`` or ``"WT|lfc"``).
        y_col : str
            Column name for the Y-axis.
        data_type : str
            ``"mle"`` (beta scores, optional TMM normalisation) or
            ``"rra"`` (log fold-change table).
        gene_col : str
            Column holding gene identifiers (default ``"id"``).
        comparison_name : str | None
            Label used in output file names (defaults to
            ``"<x_col>_vs_<y_col>"``).
        normalize : bool
            Apply ``MAGeCKFlute::NormalizeBeta`` for MLE data (default
            ``True``).
        top : int
            Number of top-label genes per selected group (default ``10``).
        groups : Sequence[str] | None
            Quadrant(s) to select hits from, e.g.
            ``["bottomleft", "topleft"]``.  Defaults to
            ``["bottomleft"]``.
        select : str | None
            Diagonal-based selection: ``"positive"``, ``"negative"``,
            ``"both"``, ``"none"``, or ``None`` (no filter).
        neg_effect_cutoff, pos_effect_cutoff : float
            Axis thresholds that define the nine quadrants.
        delta_cutoff_k : float
            Number of MADs for the diagonal cutoff (default ``2``).
        filter_fdr_x, filter_fdr_y : bool
            Require FDR significance in the respective axis column.
        filter_groups : bool
            Apply the *groups* quadrant filter when selecting hits.
        fdr_cutoff : float
            FDR significance threshold (default ``0.05``).
        toplabels : bool
            Show gene labels for selected hits (default ``True``).
        label_selected_only : bool
            Only label the selected hits, not all genes (default ``False``).
        xlab, ylab : str | None
            Axis labels (defaults to the column names).
        jpeg_width, jpeg_height : int
            Plot dimensions in cm.
        auto_cut_diag, auto_cut_x, auto_cut_y : float
            Auto-cut parameters forwarded to ``ScatterView``.
        x_cut, y_cut : tuple[float, float] | None
            Manual axis cut-off pair for RRA mode.
        tag, outdir
            Standard pipeline parameters.

        Returns
        -------
        Element
            Artefacts: ``plot`` (JPEG), ``data`` (TSV), ``qc_plot`` (JPEG),
            ``hits`` (TSV).
        """
        comparison_name = comparison_name or f"{x_col}_vs_{y_col}"
        groups = list(groups) if groups is not None else ["bottomleft"]

        # resolve input TSV from artefacts
        input_tsv: Path = table.artifacts.get("tsv") or table.artifacts.get(
            list(table.artifacts.keys())[0]
        )

        tag = from_prior(
            table.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.MAGECK,
            state=State.STAT,
            ext="jpeg",
        )

        outdir = Path(outdir or input_tsv.parent / "mageck_scatter")
        stem = f"{table.name}.scatter.{comparison_name}"
        params_path = outdir / f"{stem}.params.json"
        plot_path = outdir / f"{stem}.jpeg"
        data_path = outdir / f"{stem}_data.tsv"
        qc_plot_path = outdir / f"{stem}_lmfit.jpeg"
        hits_path = outdir / f"{stem}_hits_selected.tsv"

        payload = {
            "mode": "scatter_view",
            "input_file": input_tsv,
            "data_type": data_type,
            "output_dir": outdir,
            "filebase_name": stem,
            "x_col": x_col,
            "y_col": y_col,
            "gene_col": gene_col,
            "normalize": normalize,
            "top": top,
            "groups": groups,
            "select": select,
            "neg_effect_cutoff": neg_effect_cutoff,
            "pos_effect_cutoff": pos_effect_cutoff,
            "delta_cutoff_k": delta_cutoff_k,
            "filter_fdr_x": filter_fdr_x,
            "filter_fdr_y": filter_fdr_y,
            "filter_groups": filter_groups,
            "fdr_cutoff": fdr_cutoff,
            "toplabels": toplabels,
            "label_selected_only": label_selected_only,
            "xlab": xlab,
            "ylab": ylab,
            "jpeg_width": jpeg_width,
            "jpeg_height": jpeg_height,
            "auto_cut_diag": auto_cut_diag,
            "auto_cut_x": auto_cut_x,
            "auto_cut_y": auto_cut_y,
            "x_cut": list(x_cut) if x_cut else None,
            "y_cut": list(y_cut) if y_cut else None,
        }
        self.write_params_json(params_path, payload)

        runner = self.run_r_script(
            script_name="mageck.R",
            params_path=params_path,
            output_files=[plot_path, data_path, qc_plot_path, hits_path],
        )

        return Element(
            key=f"{table.name}.mageck.scatter.{comparison_name}.{self.version_name}",
            run=runner,
            tag=tag,
            determinants=(
                comparison_name,
                x_col,
                y_col,
                data_type,
                normalize,
                str(sorted(groups)),
                select or "",
                filter_fdr_x,
                filter_fdr_y,
                self.version_name,
            ),
            inputs=(input_tsv,),
            artifacts={
                "plot": plot_path,
                "data": data_path,
                "qc_plot": qc_plot_path,
                "hits": hits_path,
            },
            pres=(table,),
            name=f"{table.name}.mageck_scatter_{comparison_name}",
        )

    @element
    def rra(
        self,
        table: Element,
        *,
        x_col: str,
        y_col: str,
        data_type: str = "mle",
        gene_col: str = "id",
        comparison_name: str | None = None,
        normalize: bool = True,
        top: int = 10,
        groups: Sequence[str] | None = None,
        select: str | None = None,
        neg_effect_cutoff: float = -0.4,
        pos_effect_cutoff: float = 0.4,
        delta_cutoff_k: float = 2.0,
        filter_fdr_x: bool = False,
        filter_fdr_y: bool = False,
        filter_groups: bool = True,
        fdr_cutoff: float = 0.05,
        toplabels: bool = True,
        label_selected_only: bool = False,
        xlab: str | None = None,
        ylab: str | None = None,
        jpeg_width: int = 20,
        jpeg_height: int = 15,
        auto_cut_diag: float = 2.0,
        auto_cut_x: float = 2.0,
        auto_cut_y: float = 2.0,
        x_cut: tuple[float, float] | None = None,
        y_cut: tuple[float, float] | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
    ) -> Element:
        """Scatter-view hit calling via MAGeCKFlute.

        Reads a MAGeCK MLE or RRA result TSV, produces a ``ScatterView``
        JPEG, a QC linear-fit JPEG, a data TSV with enrichment annotations,
        and a TSV of selected hit genes.

        Parameters
        ----------
        table : Element
            Input element.  Must expose a ``tsv`` artefact pointing to the
            MAGeCK result table.
        x_col : str
            Column name for the X-axis (e.g. ``"WT|beta"`` or ``"WT|lfc"``).
        y_col : str
            Column name for the Y-axis.
        data_type : str
            ``"mle"`` (beta scores, optional TMM normalisation) or
            ``"rra"`` (log fold-change table).
        gene_col : str
            Column holding gene identifiers (default ``"id"``).
        comparison_name : str | None
            Label used in output file names (defaults to
            ``"<x_col>_vs_<y_col>"``).
        normalize : bool
            Apply ``MAGeCKFlute::NormalizeBeta`` for MLE data (default
            ``True``).
        top : int
            Number of top-label genes per selected group (default ``10``).
        groups : Sequence[str] | None
            Quadrant(s) to select hits from, e.g.
            ``["bottomleft", "topleft"]``.  Defaults to
            ``["bottomleft"]``.
        select : str | None
            Diagonal-based selection: ``"positive"``, ``"negative"``,
            ``"both"``, ``"none"``, or ``None`` (no filter).
        neg_effect_cutoff, pos_effect_cutoff : float
            Axis thresholds that define the nine quadrants.
        delta_cutoff_k : float
            Number of MADs for the diagonal cutoff (default ``2``).
        filter_fdr_x, filter_fdr_y : bool
            Require FDR significance in the respective axis column.
        filter_groups : bool
            Apply the *groups* quadrant filter when selecting hits.
        fdr_cutoff : float
            FDR significance threshold (default ``0.05``).
        toplabels : bool
            Show gene labels for selected hits (default ``True``).
        label_selected_only : bool
            Only label the selected hits, not all genes (default ``False``).
        xlab, ylab : str | None
            Axis labels (defaults to the column names).
        jpeg_width, jpeg_height : int
            Plot dimensions in cm.
        auto_cut_diag, auto_cut_x, auto_cut_y : float
            Auto-cut parameters forwarded to ``ScatterView``.
        x_cut, y_cut : tuple[float, float] | None
            Manual axis cut-off pair for RRA mode.
        tag, outdir
            Standard pipeline parameters.

        Returns
        -------
        Element
            Artefacts: ``plot`` (JPEG), ``data`` (TSV), ``qc_plot`` (JPEG),
            ``hits`` (TSV).
        """
        comparison_name = comparison_name or f"{x_col}_vs_{y_col}"
        groups = list(groups) if groups is not None else ["bottomleft"]

        # resolve input TSV from artefacts
        input_tsv: Path = table.artifacts.get("tsv") or table.artifacts.get(
            list(table.artifacts.keys())[0]
        )

        tag = from_prior(
            table.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.MAGECK,
            state=State.STAT,
            ext="jpeg",
        )

        outdir = Path(outdir or input_tsv.parent / "mageck_scatter")
        stem = f"{table.name}.scatter.{comparison_name}"
        params_path = outdir / f"{stem}.params.json"
        plot_path = outdir / f"{stem}.jpeg"
        data_path = outdir / f"{stem}_data.tsv"
        qc_plot_path = outdir / f"{stem}_lmfit.jpeg"
        hits_path = outdir / f"{stem}_hits_selected.tsv"

        payload = {
            "mode": "scatter_view",
            "input_file": input_tsv,
            "data_type": data_type,
            "output_dir": outdir,
            "filebase_name": stem,
            "x_col": x_col,
            "y_col": y_col,
            "gene_col": gene_col,
            "normalize": normalize,
            "top": top,
            "groups": groups,
            "select": select,
            "neg_effect_cutoff": neg_effect_cutoff,
            "pos_effect_cutoff": pos_effect_cutoff,
            "delta_cutoff_k": delta_cutoff_k,
            "filter_fdr_x": filter_fdr_x,
            "filter_fdr_y": filter_fdr_y,
            "filter_groups": filter_groups,
            "fdr_cutoff": fdr_cutoff,
            "toplabels": toplabels,
            "label_selected_only": label_selected_only,
            "xlab": xlab,
            "ylab": ylab,
            "jpeg_width": jpeg_width,
            "jpeg_height": jpeg_height,
            "auto_cut_diag": auto_cut_diag,
            "auto_cut_x": auto_cut_x,
            "auto_cut_y": auto_cut_y,
            "x_cut": list(x_cut) if x_cut else None,
            "y_cut": list(y_cut) if y_cut else None,
        }
        self.write_params_json(params_path, payload)

        runner = self.mageck_test(
            script_name="mageck.R",
            params_path=params_path,
            output_files=[plot_path, data_path, qc_plot_path, hits_path],
        )

        return Element(
            key=f"{table.name}.mageck.scatter.{comparison_name}.{self.version_name}",
            run=runner,
            tag=tag,
            determinants=(
                comparison_name,
                x_col,
                y_col,
                data_type,
                normalize,
                str(sorted(groups)),
                select or "",
                filter_fdr_x,
                filter_fdr_y,
                self.version_name,
            ),
            inputs=(input_tsv,),
            artifacts={
                "plot": plot_path,
                "data": data_path,
                "qc_plot": qc_plot_path,
                "hits": hits_path,
            },
            pres=(table,),
            name=f"{table.name}.mageck_scatter_{comparison_name}",
        )
