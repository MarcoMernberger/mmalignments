"""edgeR differential expression — pipeline Element interface.

Each public method decorated with ``@element`` accepts a :class:`TableElement`
(raw count matrix) and returns a :class:`TableElement` with per-gene statistics.

The actual computation is delegated to ``Rscript edger.R`` via
:meth:`RScriptExternal.run_r_script`.

Available analyses
------------------
* :meth:`unpaired` — exactTest two-group comparison with TMM normalisation.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Mapping, Sequence

from mmalignments.models.elements import TableElement, element
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

_R_SCRIPT_DIR = Path(__file__).parent.parent.parent / "r"

_DEFAULT_COLUMN_MAP_UNPAIRED: dict[str, str] = {
    "logFC":   "log2FC ({comparison})",
    "PValue":  "p ({comparison})",
    "FDR":     "FDR ({comparison})",
    "logCPM":  "logCPM ({comparison})",
}


class EdgeR(RScriptExternal):
    """Pipeline interface for edgeR differential expression analysis.

    Examples
    --------
    >>> edger = EdgeR()
    >>> result = edger.unpaired(
    ...     counts_element,
    ...     condition_a="KO",
    ...     condition_b="WT",
    ...     columns_a=["KO_1", "KO_2"],
    ...     columns_b=["WT_1", "WT_2"],
    ...     comparison_name="KO_vs_WT",
    ... )
    """

    def __init__(
        self,
        r_script_dir: Path | None = None,
        version: str | None = None,
    ) -> None:
        super().__init__(
            name="EdgeR",
            r_script_dir=r_script_dir or _R_SCRIPT_DIR,
            version=version,
            source="https://bioconductor.org/packages/edgeR",
        )

    @element
    def unpaired(
        self,
        counts: TableElement,
        *,
        condition_a: str,
        condition_b: str,
        columns_a: Sequence[str],
        columns_b: Sequence[str],
        comparison_name: str | None = None,
        library_sizes: Sequence[float] | None = None,
        manual_dispersion_value: float = 0.4,
        column_map: Mapping[str, str] | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
    ) -> TableElement:
        """Two-group exactTest with TMM normalisation via edgeR.

        Parameters
        ----------
        counts : TableElement
            Raw integer count matrix (genes × samples) as Parquet.
        condition_a : str
            Treatment condition label (numerator of fold change).
        condition_b : str
            Reference condition label (denominator).
        columns_a : Sequence[str]
            Sample column names belonging to *condition_a*.
        columns_b : Sequence[str]
            Sample column names belonging to *condition_b*.
        comparison_name : str | None
            Label for file names and column names.
            Defaults to ``"<condition_a>_vs_<condition_b>"``.
        library_sizes : Sequence[float] | None
            Override library sizes (e.g. from spike-ins).  If ``None``
            column sums are used.
        manual_dispersion_value : float
            Fallback common dispersion for single-replicate comparisons
            (default ``0.4``).
        column_map : Mapping[str, str] | None
            Custom renaming of edgeR output columns.
        tag, outdir
            Standard pipeline parameters.

        Returns
        -------
        TableElement
            Artefacts: ``tsv`` and ``parquet``.
        """
        comparison_name = comparison_name or f"{condition_a}_vs_{condition_b}"

        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.DIFF,
            method=Method.EDGER,
            state=State.DIFF,
            ext="tsv",
        )

        outdir = Path(outdir or counts.tsv.parent / "edger")
        stem = f"{counts.name}.edger.{comparison_name}"
        tsv_path     = outdir / f"{stem}.tsv"
        parquet_path = outdir / f"{stem}.parquet"
        params_path  = outdir / f"{stem}.params.json"

        col_map = dict(column_map) if column_map else {
            k: v.format(comparison=comparison_name)
            for k, v in _DEFAULT_COLUMN_MAP_UNPAIRED.items()
        }

        payload = {
            "counts_parquet":           counts.parquet,
            "condition_a":              condition_a,
            "condition_b":              condition_b,
            "columns_a":                list(columns_a),
            "columns_b":                list(columns_b),
            "comparison_name":          comparison_name,
            "library_sizes":            list(library_sizes) if library_sizes else None,
            "manual_dispersion_value":  manual_dispersion_value,
            "column_map":               col_map,
            "output_tsv":               tsv_path,
            "output_parquet":           parquet_path,
        }
        self.write_params_json(params_path, payload)

        runner = self.run_r_script(
            script_name="edger.R",
            params_path=params_path,
            output_files=[tsv_path, parquet_path],
        )

        return TableElement(
            key=f"{counts.name}.edger.unpaired.{comparison_name}.{self.version_name}",
            run=runner,
            tag=tag,
            tsv=tsv_path,
            parquet=parquet_path,
            determinants=(
                comparison_name,
                condition_a,
                condition_b,
                manual_dispersion_value,
                self.version_name,
            ),
            inputs=(counts.parquet,),
            pres=(counts,),
            name=f"{counts.name}.edger_{comparison_name}",
        )

    # ------------------------------------------------------------------
    # TMM normalisation
    # ------------------------------------------------------------------

    @element
    def tmm(
        self,
        counts: TableElement,
        *,
        sample_columns: Sequence[str] | None = None,
        library_sizes: Sequence[float] | None = None,
        log: bool = True,
        prior_count: float = 2.0,
        comparison_name: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
    ) -> TableElement:
        """TMM normalisation → log2-CPM matrix via edgeR.

        Applies ``calcNormFactors(method=\"TMM\")`` followed by
        ``cpm(log=TRUE)`` and writes the resulting genes × samples matrix.

        Parameters
        ----------
        counts : TableElement
            Raw integer count matrix (genes × samples) as Parquet.
        sample_columns : Sequence[str] | None
            Columns to include.  Defaults to all columns.
        library_sizes : Sequence[float] | None
            Override library sizes (e.g. spike-in based).  If ``None``
            column sums are used.
        log : bool
            Return log2-CPM (default ``True``).  Set ``False`` for raw CPM.
        prior_count : float
            Prior count added to each observation before log transformation
            (default ``2``, avoids log(0)).
        comparison_name : str | None
            Label for file names (default ``"tmm"``).
        tag, outdir
            Standard pipeline parameters.

        Returns
        -------
        TableElement
            Artefacts: ``tsv`` and ``parquet`` — genes × samples matrix of
            TMM-normalised (log2-)CPM values.
        """
        comparison_name = comparison_name or "tmm"

        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.EDGER,
            state=State.NORMAL,
            ext="tsv",
        )

        outdir = Path(outdir or counts.tsv.parent / "edger")
        stem = f"{counts.name}.edger.{comparison_name}"
        tsv_path     = outdir / f"{stem}.tsv"
        parquet_path = outdir / f"{stem}.parquet"
        params_path  = outdir / f"{stem}.params.json"

        payload: dict = {
            "mode":           "tmm",
            "counts_parquet": counts.parquet,
            "log":            log,
            "prior_count":    prior_count,
            "output_tsv":     tsv_path,
            "output_parquet": parquet_path,
        }
        if sample_columns is not None:
            payload["sample_columns"] = list(sample_columns)
        if library_sizes is not None:
            payload["library_sizes"] = list(library_sizes)

        self.write_params_json(params_path, payload)

        runner = self.run_r_script(
            script_name="edger.R",
            params_path=params_path,
            output_files=[tsv_path, parquet_path],
        )

        # column_roles: TMM columns + propagated annotations from input
        cols = list(sample_columns) if sample_columns is not None else []
        tmm_role = "log2cpm" if log else "cpm"
        column_roles: dict[str, str] = {c: tmm_role for c in cols}
        column_roles.update(counts.propagate())

        return TableElement(
            key=f"{counts.name}.edger.tmm.{comparison_name}.{self.version_name}",
            run=runner,
            tag=tag,
            tsv=tsv_path,
            parquet=parquet_path,
            column_roles=column_roles or None,
            determinants=(
                comparison_name,
                log,
                prior_count,
                self.version_name,
            ),
            inputs=(counts.parquet,),
            pres=(counts,),
            name=f"{counts.name}.edger_{comparison_name}",
        )
