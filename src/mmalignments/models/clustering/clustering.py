"""Unsupervised clustering, PCA, and visualisation — pipeline Element interface.

Each public method decorated with ``@element`` accepts a :class:`TableElement`
(count / expression matrix) and returns an :class:`Element` whose artefacts are
plot files (HTML and/or PDF) and, where applicable, a Parquet of computed values.

The actual computation is delegated to ``Rscript clustering.R`` via
:meth:`RScriptExternal.run_r_script`.  All parameters are passed as a JSON file so
no rpy2 dependency is needed at the Python level.

Available analyses
------------------
* :meth:`pca`     — PCA scatter plot  (interactive HTML + static PDF + coords Parquet).
                    Interactive HTML: label-toggle button and font-size slider for
                    label spacing.  Static PDF: ggplot2 + ggrepel with configurable
                    ``box.padding``.
* :meth:`scree`   — Scree and elbow plots (PDF) that show variance explained per PC.
* :meth:`heatmap` — Hierarchical-clustering heatmap of the most-variable genes (PDF).
* :meth:`volcano` — Volcano plot from differential expression results
                    (interactive HTML + static PDF).
"""

from __future__ import annotations

import logging
from collections.abc import Collection
from pathlib import Path
from typing import Sequence

from mmalignments.models.elements import Element, TableElement, element
from mmalignments.models.tags import (
    ElementTag,
    Method,
    PartialElementTag,
    Stage,
    State,
    from_prior,
)
from mmalignments.r.r_integration import RScriptExternal

logger = logging.getLogger(__name__)

# R scripts live in  mmalignments/r/  alongside deseq2.R / edger.R
_R_SCRIPT_DIR = Path(__file__).parent.parent.parent / "r"


class Clustering(RScriptExternal):
    """Pipeline interface for unsupervised clustering and visualisation.

    Generates PCA scatter plots, scree/elbow plots, hierarchical-clustering
    heatmaps, and volcano plots.  All heavy lifting is done in R via a single
    ``clustering.R`` script dispatched by a ``mode`` field in the params JSON.

    Examples
    --------
    >>> cl = Clustering()
    >>> pca_elem = cl.pca(
    ...     counts_element,
    ...     metadata_csv=Path("metadata.csv"),
    ...     color_by="condition",
    ...     comparison_name="all_samples",
    ... )
    >>> scree_elem = cl.scree(counts_element)
    >>> heatmap_elem = cl.heatmap(counts_element, metadata_csv=Path("metadata.csv"),
    ...                           color_by="condition", top_n_genes=500)
    >>> volcano_elem = cl.volcano(
    ...     deseq2_result_element,
    ...     logfc_col="log2FC (KO_vs_WT)",
    ...     fdr_col="FDR (KO_vs_WT)",
    ...     comparison_name="KO_vs_WT",
    ... )
    """

    def __init__(
        self,
        r_script_dir: Path | None = None,
        version: str | None = None,
    ) -> None:
        super().__init__(
            name="Clustering",
            r_script_dir=r_script_dir or _R_SCRIPT_DIR,
            version=version,
            source="https://www.r-project.org",
        )

    # ------------------------------------------------------------------
    # PCA scatter plot
    # ------------------------------------------------------------------

    @element
    def pca(
        self,
        counts: TableElement,
        *,
        metadata_csv: Path | str | None = None,
        color_by: str | None = None,
        shape_by: str | None = None,
        show_labels: bool = True,
        label_spacing: float = 1.0,
        pcs: Sequence[tuple[int, int]] | None = None,
        n_components: int = 50,
        scale: bool = True,
        center: bool = True,
        comparison_name: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
    ) -> Element:
        """PCA scatter plot of samples.

        Runs PCA on the transposed count matrix (samples × genes) and produces:

        * **HTML** — interactive Plotly scatter with a *Show / Hide Labels* toggle
          button and a font-size slider that adjusts label density.
        * **PDF** — static ggplot2 + ggrepel scatter; ``label_spacing`` controls
          ``ggrepel::box.padding`` (larger = labels pushed further apart).
        * **Parquet** — PCA coordinates (samples × PCs) for downstream use.

        Multiple PC pairs can be requested via *pcs*; they appear as separate
        panels in the PDF and as separate tabs / a subplot in the HTML.

        Parameters
        ----------
        counts : TableElement
            Raw or normalised count matrix (genes × samples).
        metadata_csv : Path | str | None
            Optional CSV file — one row per sample, first column must match the
            sample column names in the count matrix.  Extra columns are available
            for *color_by* / *shape_by*.
        color_by : str | None
            Metadata column mapped to point colour.
        shape_by : str | None
            Metadata column mapped to point shape.
        show_labels : bool
            Whether to show sample labels when the plot first opens (default
            ``True``).  The HTML slider allows toggling this interactively.
        label_spacing : float
            ``ggrepel::box.padding`` for the static PDF.  Larger values push
            labels further apart (default ``1.0``).
        pcs : Sequence[tuple[int, int]] | None
            PC pairs to plot, e.g. ``[(1, 2), (1, 3), (2, 3)]``.
            Defaults to ``[(1, 2)]``.
        n_components : int
            Number of principal components to compute (default ``50``).
        scale : bool
            Scale each gene to unit variance before PCA (default ``True``).
        center : bool
            Center each gene before PCA (default ``True``).
        comparison_name : str | None
            Short label used in file names and plot titles (default ``"pca"``).
        tag : PartialElementTag | ElementTag | None
            Override for the output element tag.
        outdir : Path | str | None
            Output directory.  Defaults to ``<counts_dir>/clustering``.

        Returns
        -------
        Element
            Artefacts: ``html`` (interactive), ``pdf`` (static), ``parquet``
            (PCA coordinates).
        """
        comparison_name = comparison_name or "pca"
        pcs = list(pcs) if pcs is not None else [(1, 2)]

        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.CLUSTERING,
            state=State.PCA,
            ext="html",
        )

        outdir = Path(outdir or counts.tsv.parent / "clustering")
        stem = f"{counts.name}.pca.{comparison_name}"
        html_path    = outdir / f"{stem}.html"
        pdf_path     = outdir / f"{stem}.pdf"
        parquet_path = outdir / f"{stem}.parquet"
        params_path  = outdir / f"{stem}.params.json"

        payload = {
            "mode":            "pca",
            "counts_parquet":  counts.parquet,
            "metadata_csv":    metadata_csv,
            "color_by":        color_by,
            "shape_by":        shape_by,
            "show_labels":     show_labels,
            "label_spacing":   label_spacing,
            "pcs":             [list(pair) for pair in pcs],
            "n_components":    n_components,
            "scale":           scale,
            "center":          center,
            "comparison_name": comparison_name,
            "output_html":     html_path,
            "output_pdf":      pdf_path,
            "output_parquet":  parquet_path,
        }
        self.write_params_json(params_path, payload)

        runner = self.run_r_script(
            script_name="clustering.R",
            params_path=params_path,
            output_files=[html_path, pdf_path, parquet_path],
        )

        return Element(
            key=f"{counts.name}.pca.{comparison_name}.{self.version_name}",
            run=runner,
            tag=tag,
            determinants=(
                comparison_name,
                color_by or "",
                shape_by or "",
                show_labels,
                label_spacing,
                n_components,
                scale,
                center,
                self.version_name,
            ),
            inputs=(counts.parquet,),
            artifacts={"html": html_path, "pdf": pdf_path, "parquet": parquet_path},
            pres=(counts,),
            name=f"{counts.name}.pca_{comparison_name}",
        )

    # ------------------------------------------------------------------
    # Scree / elbow plot
    # ------------------------------------------------------------------

    @element
    def scree(
        self,
        counts: TableElement,
        *,
        n_components: int = 30,
        scale: bool = True,
        center: bool = True,
        comparison_name: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
    ) -> Element:
        """Scree and elbow plots of PCA variance explained.

        Runs PCA on the transposed count matrix and produces a two-panel PDF:

        * **Scree plot** — bar chart of variance explained (%) per PC with a
          vertical dashed line marking the estimated elbow point.
        * **Cumulative variance** — line chart with the cumulative variance
          explained curve; a horizontal dotted line marks the 80 % threshold.

        The elbow is estimated with a simple kneedle approximation: the PC at
        the maximum distance from the straight line between the first and last
        point on the scree curve.

        Parameters
        ----------
        counts : TableElement
            Count / expression matrix element.
        n_components : int
            Number of PCs to display (default ``30``).
        scale, center : bool
            Passed to ``prcomp`` (default both ``True``).
        comparison_name : str | None
            Label used in file names and titles (default ``"scree"``).
        tag, outdir
            Standard pipeline parameters.

        Returns
        -------
        Element
            Artefacts: ``pdf``.
        """
        comparison_name = comparison_name or "scree"

        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.CLUSTERING,
            state=State.PCA,
            ext="pdf",
        )

        outdir = Path(outdir or counts.tsv.parent / "clustering")
        stem = f"{counts.name}.scree.{comparison_name}"
        pdf_path    = outdir / f"{stem}.pdf"
        params_path = outdir / f"{stem}.params.json"

        payload = {
            "mode":            "scree",
            "counts_parquet":  counts.parquet,
            "n_components":    n_components,
            "scale":           scale,
            "center":          center,
            "comparison_name": comparison_name,
            "output_pdf":      pdf_path,
        }
        self.write_params_json(params_path, payload)

        runner = self.run_r_script(
            script_name="clustering.R",
            params_path=params_path,
            output_files=[pdf_path],
        )

        return Element(
            key=f"{counts.name}.scree.{comparison_name}.{self.version_name}",
            run=runner,
            tag=tag,
            determinants=(
                comparison_name,
                n_components,
                scale,
                center,
                self.version_name,
            ),
            inputs=(counts.parquet,),
            artifacts={"pdf": pdf_path},
            pres=(counts,),
            name=f"{counts.name}.scree_{comparison_name}",
        )

    # ------------------------------------------------------------------
    # Hierarchical-clustering heatmap
    # ------------------------------------------------------------------

    @element
    def heatmap(
        self,
        counts: TableElement,
        *,
        metadata_csv: Path | str | None = None,
        color_by: str | Collection[str] | None = None,
        top_n_genes: int = 500,
        clustering_distance_rows: str = "euclidean",
        clustering_distance_cols: str = "euclidean",
        clustering_method: str = "complete",
        scale_rows: bool = True,
        comparison_name: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
    ) -> Element:
        """Hierarchical-clustering heatmap of the most-variable genes.

        Selects the *top_n_genes* most variable genes by inter-quartile range,
        optionally z-score scales rows, and produces a clustered heatmap (PDF)
        using ``pheatmap``.  Sample annotation colour bars are added when
        *metadata_csv* and *color_by* are provided.

        Gene row names are shown only when ``top_n_genes ≤ 50``; for larger
        matrices they are suppressed to keep the plot readable.

        Parameters
        ----------
        counts : TableElement
            Count / expression matrix element (genes × samples).
        metadata_csv : Path | str | None
            Optional CSV with sample metadata — first column must match sample
            columns in the count matrix.
        color_by : str | Collection[str] | None
            Column(s) in the metadata to use as colour-bar annotation.  A single
            string or a list of strings are both accepted.
        top_n_genes : int
            Number of most-variable genes to include (default ``500``).
        clustering_distance_rows : str
            Row-clustering distance metric passed to ``pheatmap``
            (default ``"euclidean"``).
        clustering_distance_cols : str
            Column-clustering distance metric (default ``"euclidean"``).
        clustering_method : str
            Linkage method passed to ``pheatmap`` (default ``"complete"``).
        scale_rows : bool
            Z-score scale rows before plotting (default ``True``).
        comparison_name : str | None
            Label used in file names and titles (default ``"heatmap"``).
        tag, outdir
            Standard pipeline parameters.

        Returns
        -------
        Element
            Artefacts: ``pdf``.
        """
        comparison_name = comparison_name or "heatmap"
        color_by_list = (
            [color_by]
            if isinstance(color_by, str)
            else list(color_by)
            if color_by is not None
            else None
        )

        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.CLUSTERING,
            state=State.CLUSTER,
            ext="pdf",
        )

        outdir = Path(outdir or counts.tsv.parent / "clustering")
        stem = f"{counts.name}.heatmap.{comparison_name}"
        pdf_path    = outdir / f"{stem}.pdf"
        params_path = outdir / f"{stem}.params.json"

        payload = {
            "mode":                       "heatmap",
            "counts_parquet":             counts.parquet,
            "metadata_csv":               metadata_csv,
            "color_by":                   color_by_list,
            "top_n_genes":                top_n_genes,
            "clustering_distance_rows":   clustering_distance_rows,
            "clustering_distance_cols":   clustering_distance_cols,
            "clustering_method":          clustering_method,
            "scale_rows":                 scale_rows,
            "comparison_name":            comparison_name,
            "output_pdf":                 pdf_path,
        }
        self.write_params_json(params_path, payload)

        runner = self.run_r_script(
            script_name="clustering.R",
            params_path=params_path,
            output_files=[pdf_path],
        )

        return Element(
            key=f"{counts.name}.heatmap.{comparison_name}.{self.version_name}",
            run=runner,
            tag=tag,
            determinants=(
                comparison_name,
                top_n_genes,
                clustering_distance_rows,
                clustering_distance_cols,
                clustering_method,
                scale_rows,
                self.version_name,
            ),
            inputs=(counts.parquet,),
            artifacts={"pdf": pdf_path},
            pres=(counts,),
            name=f"{counts.name}.heatmap_{comparison_name}",
        )

    # ------------------------------------------------------------------
    # Volcano plot
    # ------------------------------------------------------------------

    @element
    def volcano(
        self,
        results: TableElement,
        *,
        logfc_col: str,
        fdr_col: str,
        p_col: str | None = None,
        logfc_threshold: float = 1.0,
        fdr_threshold: float = 0.05,
        label_top_n: int = 20,
        label_col: str | None = None,
        label_spacing: float = 1.0,
        show_labels: bool = True,
        comparison_name: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
    ) -> Element:
        """Volcano plot from differential expression results.

        Plots −log₁₀(FDR) vs log₂FC and colours points by significance
        category (up / down / n.s.) based on *fdr_threshold* and
        *logfc_threshold*.  The top *label_top_n* genes (ranked by
        FDR × |log₂FC|) are labelled.

        Produces:

        * **HTML** — interactive Plotly scatter with rich hover tooltips, a
          *Show / Hide Labels* toggle button, and a font-size slider for
          adjusting label density.
        * **PDF** — static ggplot2 + ggrepel plot; ``label_spacing`` controls
          ``ggrepel::box.padding``.

        Parameters
        ----------
        results : TableElement
            Differential expression result element (e.g. from :class:`DESeq2`
            or :class:`EdgeR`).  Must contain *logfc_col* and *fdr_col*.
        logfc_col : str
            Column name for log₂ fold-change values.
        fdr_col : str
            Column name for adjusted p-values / FDR values.  Used on the y-axis
            unless *p_col* is given.
        p_col : str | None
            Column name for raw p-values to display on the y-axis instead of
            FDR.  If ``None``, *fdr_col* is used for both colouring and y-axis.
        logfc_threshold : float
            |log₂FC| threshold for colouring points (default ``1.0``).
        fdr_threshold : float
            FDR threshold for colouring points (default ``0.05``).
        label_top_n : int
            Number of top-genes to label (default ``20``).
        label_col : str | None
            Column to use as gene labels.  Defaults to the row-index of the
            result table.
        label_spacing : float
            ``ggrepel::box.padding`` value for the static PDF (default ``1.0``).
        show_labels : bool
            Whether labels are shown when the HTML first opens (default ``True``).
        comparison_name : str | None
            Label used in file names and plot titles (default ``"volcano"``).
        tag, outdir
            Standard pipeline parameters.

        Returns
        -------
        Element
            Artefacts: ``html`` (interactive), ``pdf`` (static).
        """
        comparison_name = comparison_name or "volcano"

        tag = from_prior(
            results.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.CLUSTERING,
            state=State.DIFF,
            ext="html",
        )

        outdir = Path(outdir or results.tsv.parent / "clustering")
        stem = f"{results.name}.volcano.{comparison_name}"
        html_path   = outdir / f"{stem}.html"
        pdf_path    = outdir / f"{stem}.pdf"
        params_path = outdir / f"{stem}.params.json"

        payload = {
            "mode":             "volcano",
            "results_parquet":  results.parquet,
            "logfc_col":        logfc_col,
            "fdr_col":          fdr_col,
            "p_col":            p_col,
            "logfc_threshold":  logfc_threshold,
            "fdr_threshold":    fdr_threshold,
            "label_top_n":      label_top_n,
            "label_col":        label_col,
            "label_spacing":    label_spacing,
            "show_labels":      show_labels,
            "comparison_name":  comparison_name,
            "output_html":      html_path,
            "output_pdf":       pdf_path,
        }
        self.write_params_json(params_path, payload)

        runner = self.run_r_script(
            script_name="clustering.R",
            params_path=params_path,
            output_files=[html_path, pdf_path],
        )

        return Element(
            key=f"{results.name}.volcano.{comparison_name}.{self.version_name}",
            run=runner,
            tag=tag,
            determinants=(
                comparison_name,
                logfc_col,
                fdr_col,
                logfc_threshold,
                fdr_threshold,
                label_top_n,
                label_spacing,
                self.version_name,
            ),
            inputs=(results.parquet,),
            artifacts={"html": html_path, "pdf": pdf_path},
            pres=(results,),
            name=f"{results.name}.volcano_{comparison_name}",
        )
