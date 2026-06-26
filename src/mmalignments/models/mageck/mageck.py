"""Wrapper for the MAGeCK CLI tool.

This module provides a pipeline-friendly interface for the Python MAGeCK CLI
(``mageck``), following the same high-level/low-level pattern used in
``mmfqcount``:

* high-level ``@element`` methods model DAG dependencies and cache signatures
* low-level ``@subroutine`` methods build the concrete CLI calls

Implemented subcommands
----------------------
* ``count``
* ``test``
* ``mle``
* ``plot``
* ``pathway``
"""

from __future__ import annotations

import logging
import re
import subprocess
from pathlib import Path
from typing import Mapping, Sequence, Callable, Iterable
from pandas import DataFrame, Series
from statsmodels.stats.multitest import multipletests
from mmalignments.models.elements import Element, TableElement, element
from mmalignments.models.parameters import (
    ParamRegistry,
    Params,
    ParamSet,
    ParamSpec,
    ToolThreadSpec,
    render_flag,
    render_value,
)
from mmalignments.models.tags import (
    ElementTag,
    Method,
    PartialElementTag,
    Stage,
    State,
    from_prior,
)

from ..externals import External, ExternalRunConfig, SubroutineIn, subroutine

logger = logging.getLogger(__name__)


def _build_param_registry() -> ParamRegistry:
    """Return MAGeCK parameter registry by subcommand."""
    count_specs: dict[str, ParamSpec] = {
        "norm_method": ParamSpec(
            "norm_method", "--norm-method", str, render=render_value
        ),
        "control_sgrna": ParamSpec(
            "control_sgrna", "--control-sgrna", str, render=render_value
        ),
        "control_gene": ParamSpec(
            "control_gene", "--control-gene", str, render=render_value
        ),
        "sample_label": ParamSpec(
            "sample_label", "--sample-label", str, render=render_value
        ),
        "unmapped_to_file": ParamSpec(
            "unmapped_to_file", "--unmapped-to-file", bool, render=render_flag
        ),
        "keep_tmp": ParamSpec("keep_tmp", "--keep-tmp", bool, render=render_flag),
        "test_run": ParamSpec("test_run", "--test-run", bool, render=render_flag),
        "count_pair": ParamSpec("count_pair", "--count-pair", bool, render=render_flag),
        "trim_5": ParamSpec("trim_5", "--trim-5", str, render=render_value),
        "sgrna_len": ParamSpec("sgrna_len", "--sgrna-len", int, render=render_value),
        "count_n": ParamSpec("count_n", "--count-n", bool, render=render_flag),
        "reverse_complement": ParamSpec(
            "reverse_complement", "--reverse-complement", bool, render=render_flag
        ),
        "pdf_report": ParamSpec("pdf_report", "--pdf-report", bool, render=render_flag),
        "day0_label": ParamSpec("day0_label", "--day0-label", str, render=render_value),
        "gmt_file": ParamSpec("gmt_file", "--gmt-file", str, render=render_value),
    }

    test_specs: dict[str, ParamSpec] = {
        "control_id": ParamSpec("control_id", "--control-id", str, render=render_value),
        "paired": ParamSpec("paired", "--paired", bool, render=render_flag),
        "norm_method": ParamSpec(
            "norm_method", "--norm-method", str, render=render_value
        ),
        "gene_test_fdr_threshold": ParamSpec(
            "gene_test_fdr_threshold",
            "--gene-test-fdr-threshold",
            float,
            render=render_value,
        ),
        "adjust_method": ParamSpec(
            "adjust_method", "--adjust-method", str, render=render_value
        ),
        "variance_estimation_samples": ParamSpec(
            "variance_estimation_samples",
            "--variance-estimation-samples",
            str,
            render=render_value,
        ),
        "sort_criteria": ParamSpec(
            "sort_criteria", "--sort-criteria", str, render=render_value
        ),
        "remove_zero": ParamSpec(
            "remove_zero", "--remove-zero", str, render=render_value
        ),
        "remove_zero_threshold": ParamSpec(
            "remove_zero_threshold",
            "--remove-zero-threshold",
            float,
            render=render_value,
        ),
        "pdf_report": ParamSpec("pdf_report", "--pdf-report", bool, render=render_flag),
        "gene_lfc_method": ParamSpec(
            "gene_lfc_method", "--gene-lfc-method", str, render=render_value
        ),
        "control_sgrna": ParamSpec(
            "control_sgrna", "--control-sgrna", str, render=render_value
        ),
        "control_gene": ParamSpec(
            "control_gene", "--control-gene", str, render=render_value
        ),
        "normcounts_to_file": ParamSpec(
            "normcounts_to_file", "--normcounts-to-file", bool, render=render_flag
        ),
        "skip_gene": ParamSpec("skip_gene", "--skip-gene", str, render=render_value),
        "keep_tmp": ParamSpec("keep_tmp", "--keep-tmp", bool, render=render_flag),
        "additional_rra_parameters": ParamSpec(
            "additional_rra_parameters",
            "--additional-rra-parameters",
            str,
            render=render_value,
        ),
        "cnv_norm": ParamSpec("cnv_norm", "--cnv-norm", str, render=render_value),
        "cell_line": ParamSpec("cell_line", "--cell-line", str, render=render_value),
        "cnv_est": ParamSpec("cnv_est", "--cnv-est", str, render=render_value),
    }

    mle_specs: dict[str, ParamSpec] = {
        "include_samples": ParamSpec(
            "include_samples", "--include-samples", str, render=render_value
        ),
        "beta_labels": ParamSpec(
            "beta_labels", "--beta-labels", str, render=render_value
        ),
        "control_sgrna": ParamSpec(
            "control_sgrna", "--control-sgrna", str, render=render_value
        ),
        "control_gene": ParamSpec(
            "control_gene", "--control-gene", str, render=render_value
        ),
        "cnv_norm": ParamSpec("cnv_norm", "--cnv-norm", str, render=render_value),
        "cell_line": ParamSpec("cell_line", "--cell-line", str, render=render_value),
        "cnv_est": ParamSpec("cnv_est", "--cnv-est", str, render=render_value),
        "debug": ParamSpec("debug", "--debug", bool, render=render_flag),
        "debug_gene": ParamSpec("debug_gene", "--debug-gene", str, render=render_value),
        "norm_method": ParamSpec(
            "norm_method", "--norm-method", str, render=render_value
        ),
        "genes_varmodeling": ParamSpec(
            "genes_varmodeling", "--genes-varmodeling", int, render=render_value
        ),
        "permutation_round": ParamSpec(
            "permutation_round", "--permutation-round", int, render=render_value
        ),
        "no_permutation_by_group": ParamSpec(
            "no_permutation_by_group",
            "--no-permutation-by-group",
            bool,
            render=render_flag,
        ),
        "max_sgrnapergene_permutation": ParamSpec(
            "max_sgrnapergene_permutation",
            "--max-sgrnapergene-permutation",
            int,
            render=render_value,
        ),
        "remove_outliers": ParamSpec(
            "remove_outliers", "--remove-outliers", bool, render=render_flag
        ),
        "threads": ParamSpec(
            "threads", "--threads", int, render=render_value, affects_output=False
        ),
        "adjust_method": ParamSpec(
            "adjust_method", "--adjust-method", str, render=render_value
        ),
        "sgrna_efficiency": ParamSpec(
            "sgrna_efficiency", "--sgrna-efficiency", str, render=render_value
        ),
        "sgrna_eff_name_column": ParamSpec(
            "sgrna_eff_name_column", "--sgrna-eff-name-column", str, render=render_value
        ),
        "sgrna_eff_score_column": ParamSpec(
            "sgrna_eff_score_column",
            "--sgrna-eff-score-column",
            str,
            render=render_value,
        ),
        "update_efficiency": ParamSpec(
            "update_efficiency", "--update-efficiency", bool, render=render_flag
        ),
        "bayes": ParamSpec("bayes", "--bayes", bool, render=render_flag),
        "ppi_prior": ParamSpec("ppi_prior", "--PPI-prior", bool, render=render_flag),
        "ppi_weighting": ParamSpec(
            "ppi_weighting", "--PPI-weighting", float, render=render_value
        ),
        "negative_control": ParamSpec(
            "negative_control", "--negative-control", str, render=render_value
        ),
        "_thread_spec": ToolThreadSpec(
            flag="threads", multi=True, fraction=1.0, max_threads=None
        ),
    }

    plot_specs: dict[str, ParamSpec] = {
        "genes": ParamSpec("genes", "--genes", str, render=render_value),
        "samples": ParamSpec("samples", "--samples", str, render=render_value),
        "norm_method": ParamSpec(
            "norm_method", "--norm-method", str, render=render_value
        ),
        "control_sgrna": ParamSpec(
            "control_sgrna", "--control-sgrna", str, render=render_value
        ),
        "control_gene": ParamSpec(
            "control_gene", "--control-gene", str, render=render_value
        ),
        "keep_tmp": ParamSpec("keep_tmp", "--keep-tmp", bool, render=render_flag),
    }

    pathway_specs: dict[str, ParamSpec] = {
        "method": ParamSpec("method", "--method", str, render=render_value),
        "single_ranking": ParamSpec(
            "single_ranking", "--single-ranking", bool, render=render_flag
        ),
        "sort_criteria": ParamSpec(
            "sort_criteria", "--sort-criteria", str, render=render_value
        ),
        "keep_tmp": ParamSpec("keep_tmp", "--keep-tmp", bool, render=render_flag),
        "ranking_column": ParamSpec(
            "ranking_column", "--ranking-column", str, render=render_value
        ),
        "ranking_column_2": ParamSpec(
            "ranking_column_2", "--ranking-column-2", str, render=render_value
        ),
        "pathway_alpha": ParamSpec(
            "pathway_alpha", "--pathway-alpha", float, render=render_value
        ),
        "permutation": ParamSpec(
            "permutation", "--permutation", int, render=render_value
        ),
    }

    return ParamRegistry(
        default=ParamSet({}, "mageck", "default"),
        by_subcommand={
            "count": ParamSet(count_specs, "mageck", "count"),
            "test": ParamSet(test_specs, "mageck", "test"),
            "mle": ParamSet(mle_specs, "mageck", "mle"),
            "plot": ParamSet(plot_specs, "mageck", "plot"),
            "pathway": ParamSet(pathway_specs, "mageck", "pathway"),
        },
    )


def _build_param_registry_as_mapping() -> dict[str, ParamSet]:
    """Return a mapping accepted by External.__init__."""
    registry = _build_param_registry()
    return {
        "default": registry.default,
        **registry.by_subcommand,
    }


class Mageck(External):
    """Wrapper for the MAGeCK CLI (``mageck`` binary)."""

    def __init__(
        self,
        name: str = "mageck",
        primary_binary: str = "mageck",
        version: str | None = None,
        source: str = "https://sourceforge.net/p/mageck/wiki/Home/",
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
    ) -> None:
        resolved_parameters: Mapping[str, ParamSet] | ParamSet = (
            parameters if parameters is not None else _build_param_registry_as_mapping()
        )
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=resolved_parameters,
        )
        self.param_registry = _build_param_registry()

    def get_version(self, fallback: str | None = None) -> str | None:
        """Detect ``mageck`` version from ``mageck --version`` output."""
        if self._version:
            return self._version
        if not self.primary_binary or not self.ensure_binary():
            return fallback
        try:
            cp = subprocess.run(
                [self.primary_binary, "--version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            if cp.returncode == 0:
                out = (cp.stdout or cp.stderr or "").strip()
                if out:
                    m = re.search(r"(\d+\.\d+(?:\.\d+)?)", out)
                    if m:
                        return m.group(1)
                    return out.splitlines()[0]
        except Exception:
            pass
        return fallback

    def default_output_dir(self, root: str) -> Path:
        """Default output root for MAGeCK runs."""
        return Path("results") / "mageck" / self.version_name / root

    def _resolve_path(
        self,
        data: Element | Path | str,
        *,
        preferred_keys: Sequence[str] = ("tsv", "count_table", "gene_summary", "path"),
    ) -> Path:
        """Resolve an input path from an Element or plain path."""
        if isinstance(data, Element):
            for key in preferred_keys:
                p = data.artifacts.get(key)
                if p:
                    return Path(p).absolute()
            if data.artifacts:
                return Path(next(iter(data.artifacts.values()))).absolute()
            raise ValueError("Input element has no artifacts to resolve a path")
        return Path(data).absolute()

    def _resolve_pres(self, data: Element | Path | str) -> tuple[Element, ...]:
        return (data,) if isinstance(data, Element) else ()

    @element
    def count(
        self,
        list_seq: Element | Path | str,
        *,
        fastq: Sequence[str | Path] | None = None,
        fastq_2: Sequence[str | Path] | None = None,
        count_table: Element | Path | str | None = None,
        sample_label: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        output_prefix: str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> TableElement:
        """Run ``mageck count`` and return the produced count table element."""
        if (fastq is None) == (count_table is None):
            raise ValueError("Provide exactly one of fastq or count_table")

        list_seq_path = self._resolve_path(list_seq, preferred_keys=("tsv", "path"))
        count_table_path = (
            self._resolve_path(
                count_table, preferred_keys=("tsv", "count_table", "path")
            )
            if count_table is not None
            else None
        )

        root = (
            Path(output_prefix).stem
            if output_prefix
            else (
                count_table_path.stem
                if count_table_path is not None
                else Path(str((fastq or ["input"])[0])).stem
            )
        )

        base_tag = (
            from_prior(
                list_seq.tag,
                tag,
                stage=Stage.QUANT,
                method=Method.MAGECK,
                state=State.COUNT,
                ext="txt",
            )
            if isinstance(list_seq, Element)
            else ElementTag(
                root=root,
                level=1,
                stage=Stage.QUANT,
                method=Method.MAGECK,
                state=State.COUNT,
                ext="txt",
            ).merge(tag)
        )

        output_dir = Path(outdir or self.default_output_dir(base_tag.root)).absolute()
        prefix = (
            output_prefix
            or f"{base_tag.root}.S{base_tag.level:02d}.{base_tag.stage}.{base_tag.method}.{base_tag.state}"
        )
        prefix_path = output_dir / prefix
        count_txt = prefix_path.with_suffix(".count.txt")
        summary_txt = prefix_path.with_suffix(".countsummary.txt")

        merged_params = (params or Params()).override(
            **({"sample_label": sample_label} if sample_label else {})
        )

        runner = self.run_count(
            list_seq=list_seq_path,
            output_prefix=prefix_path,
            fastq=fastq,
            fastq_2=fastq_2,
            count_table=count_table_path,
            params=merged_params,
            cfg=cfg,
        )

        key, name = self.build_element_name(base_tag, "count")
        determinants = self.signature_determinants(merged_params, subroutine="count")
        inputs = [list_seq_path]
        if count_table_path is not None:
            inputs.append(count_table_path)
        if fastq:
            inputs.extend(Path(f).absolute() for f in fastq)
        if fastq_2:
            inputs.extend(Path(f).absolute() for f in fastq_2)

        artifacts = {
            "tsv": count_txt,
            "count_table": count_txt,
            "count_summary": summary_txt,
            "output_prefix": prefix_path,
        }
        if merged_params.get("unmapped_to_file"):
            artifacts["unmapped"] = prefix_path.with_suffix(".unmapped.txt")

        return TableElement(
            key,
            runner,
            tag=base_tag,
            tsv=count_txt,
            artifacts=artifacts,
            determinants=determinants,
            inputs=tuple(inputs),
            pres=self._resolve_pres(list_seq)
            + ((count_table,) if isinstance(count_table, Element) else ()),
            name=name,
        )

    @subroutine
    def run_count(
        self,
        list_seq: Path | str,
        output_prefix: Path | str,
        *,
        fastq: Sequence[str | Path] | None = None,
        fastq_2: Sequence[str | Path] | None = None,
        count_table: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for ``mageck count``."""
        if (fastq is None) == (count_table is None):
            raise ValueError("Provide exactly one of fastq or count_table")

        list_seq_path = Path(list_seq).absolute()
        prefix = Path(output_prefix).absolute()

        arguments: list[str] = [
            "count",
            "--list-seq",
            str(list_seq_path),
            "--output-prefix",
            str(prefix),
        ]
        in_paths: list[Path] = [list_seq_path]

        if count_table is not None:
            count_table_path = Path(count_table).absolute()
            arguments.extend(["--count-table", str(count_table_path)])
            in_paths.append(count_table_path)
        else:
            fastq_paths = [Path(f).absolute() for f in (fastq or [])]
            arguments.append("--fastq")
            arguments.extend([str(f) for f in fastq_paths])
            in_paths.extend(fastq_paths)
            if fastq_2:
                fastq2_paths = [Path(f).absolute() for f in fastq_2]
                arguments.append("--fastq-2")
                arguments.extend([str(f) for f in fastq2_paths])
                in_paths.extend(fastq2_paths)

        cli_extras = self.to_cli(params or Params(), subroutine="count")
        arguments.extend(cli_extras)

        out_paths = [
            prefix.with_suffix(".count.txt"),
            prefix.with_suffix(".countsummary.txt"),
        ]
        if (params or Params()).get("unmapped_to_file"):
            out_paths.append(prefix.with_suffix(".unmapped.txt"))

        return (
            arguments,
            "count",
            in_paths,
            out_paths,
            None,
            None,
            None,
        )

    @element
    def rra(
        self,
        counts: Element,
        *,
        treatment_id: str,
        control_id: str | None = None,
        day0_label: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        prefix: str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Run ``mageck test`` and return gene-summary as primary table."""
        if not treatment_id:
            raise ValueError("treatment_id is required")
        if bool(control_id) == bool(day0_label):
            raise ValueError("Provide exactly one of control_id or day0_label")

        count_table_path = self._resolve_path(
            counts, preferred_keys=("count_table", "tsv", "path")
        )

        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.DIFF,
            method=Method.MAGECK,
            state=State.STAT,
            ext="txt",
        )

        output_dir = Path(outdir or self.default_output_dir(tag.root)).absolute()
        prefix = prefix or tag.default_output[:-4]  # remove .txt

        prefix_path = output_dir / prefix
        gene_summary = output_dir / f"{prefix}.gene_summary.txt"
        sgrna_summary = output_dir / f"{prefix}.sgrna_summary.txt"

        runner = self.run_test(
            count_table=count_table_path,
            treatment_id=treatment_id,
            control_id=control_id,
            day0_label=day0_label,
            output_prefix=prefix_path,
            params=params,
            cfg=cfg,
        )

        key, name = self.build_element_name(tag, "test", treatment=treatment_id)
        determinants = self.signature_determinants(params, subroutine="test")
        determinants += (f"treatment_id={treatment_id}",)
        if control_id:
            determinants += (f"control_id={control_id}",)
        if day0_label:
            determinants += (f"day0_label={day0_label}",)

        artifacts = {
            "tsv": sgrna_summary,
            "gene": gene_summary,
            "sgrna": sgrna_summary,
            "output_prefix": prefix_path,
        }
        if (params or Params()).get("normcounts_to_file"):
            artifacts["normalized"] = prefix_path.with_suffix(".normalized.txt")

        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=(count_table_path,),
            pres=self._resolve_pres(counts),
            name=name,
        )

    @subroutine
    def run_test(
        self,
        count_table: Path | str,
        treatment_id: str,
        output_prefix: Path | str,
        *,
        control_id: str | None = None,
        day0_label: str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for ``mageck test``."""
        if bool(control_id) == bool(day0_label):
            raise ValueError("Provide exactly one of control_id or day0_label")

        count_table_path = Path(count_table).absolute()
        prefix = Path(output_prefix).absolute()

        arguments = [
            "test",
            "--count-table",
            str(count_table_path),
            "--treatment-id",
            treatment_id,
            "--output-prefix",
            str(prefix),
        ]
        if control_id:
            arguments.extend(["--control-id", control_id])
        if day0_label:
            arguments.extend(["--day0-label", day0_label])

        arguments.extend(self.to_cli(params or Params(), subroutine="test"))
        out_paths = [
            prefix.with_suffix(".gene_summary.txt"),
            prefix.with_suffix(".sgrna_summary.txt"),
        ]
        if (params or Params()).get("normcounts_to_file"):
            out_paths.append(prefix.with_suffix(".normalized.txt"))

        return (
            arguments,
            "test",
            [count_table_path],
            out_paths,
            None,
            None,
            None,
        )

    @element
    def mle(
        self,
        count_table: Element | Path | str,
        *,
        design_matrix: Path | str | None = None,
        day0_label: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        output_prefix: str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> TableElement:
        """Run ``mageck mle`` and return gene-summary as primary table."""
        if bool(design_matrix) == bool(day0_label):
            raise ValueError("Provide exactly one of design_matrix or day0_label")

        count_table_path = self._resolve_path(
            count_table, preferred_keys=("count_table", "tsv", "path")
        )
        design_path = Path(design_matrix).absolute() if design_matrix else None

        base_tag = (
            from_prior(
                count_table.tag,
                tag,
                stage=Stage.DIFF,
                method=Method.MAGECK,
                state=State.MODEL,
                ext="txt",
            )
            if isinstance(count_table, Element)
            else ElementTag(
                root=count_table_path.stem,
                level=1,
                stage=Stage.DIFF,
                method=Method.MAGECK,
                state=State.MODEL,
                ext="txt",
            ).merge(tag)
        )

        output_dir = Path(outdir or self.default_output_dir(base_tag.root)).absolute()
        prefix = (
            output_prefix
            or f"{base_tag.root}.S{base_tag.level:02d}.{base_tag.stage}.{base_tag.method}.{base_tag.state}"
        )
        prefix_path = output_dir / prefix
        gene_summary = prefix_path.with_suffix(".gene_summary.txt")
        sgrna_summary = prefix_path.with_suffix(".sgrna_summary.txt")

        runner = self.run_mle(
            count_table=count_table_path,
            output_prefix=prefix_path,
            design_matrix=design_path,
            day0_label=day0_label,
            params=params,
            cfg=cfg,
        )

        key, name = self.build_element_name(base_tag, "mle")
        determinants = self.signature_determinants(params, subroutine="mle")
        if design_path:
            determinants += (f"design_matrix={design_path}",)
        if day0_label:
            determinants += (f"day0_label={day0_label}",)

        artifacts = {
            "tsv": gene_summary,
            "gene_summary": gene_summary,
            "sgrna_summary": sgrna_summary,
            "output_prefix": prefix_path,
        }

        inputs = [count_table_path]
        if design_path:
            inputs.append(design_path)

        return TableElement(
            key,
            runner,
            tag=base_tag,
            tsv=gene_summary,
            artifacts=artifacts,
            determinants=determinants,
            inputs=tuple(inputs),
            pres=self._resolve_pres(count_table),
            name=name,
        )

    @subroutine
    def run_mle(
        self,
        count_table: Path | str,
        output_prefix: Path | str,
        *,
        design_matrix: Path | str | None = None,
        day0_label: str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for ``mageck mle``."""
        if bool(design_matrix) == bool(day0_label):
            raise ValueError("Provide exactly one of design_matrix or day0_label")

        count_table_path = Path(count_table).absolute()
        prefix = Path(output_prefix).absolute()

        arguments = [
            "mle",
            "--count-table",
            str(count_table_path),
            "--output-prefix",
            str(prefix),
        ]
        in_paths: list[Path] = [count_table_path]

        if design_matrix is not None:
            design_path = Path(design_matrix).absolute()
            arguments.extend(["--design-matrix", str(design_path)])
            in_paths.append(design_path)
        if day0_label is not None:
            arguments.extend(["--day0-label", day0_label])

        arguments.extend(self.to_cli(params or Params(), subroutine="mle"))

        return (
            arguments,
            "mle",
            in_paths,
            [
                prefix.with_suffix(".gene_summary.txt"),
                prefix.with_suffix(".sgrna_summary.txt"),
            ],
            None,
            None,
            None,
        )

    @element
    def plot(
        self,
        count_table: Element | Path | str,
        gene_summary: Element | Path | str,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        output_prefix: str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Run ``mageck plot`` and return an element tracking plot artefacts."""
        count_table_path = self._resolve_path(
            count_table, preferred_keys=("count_table", "tsv", "path")
        )
        gene_summary_path = self._resolve_path(
            gene_summary, preferred_keys=("gene_summary", "tsv", "path")
        )

        base_tag = (
            from_prior(
                gene_summary.tag,
                tag,
                stage=Stage.ANALYSIS,
                method=Method.MAGECK,
                state=State.REPORT,
                ext="pdf",
            )
            if isinstance(gene_summary, Element)
            else ElementTag(
                root=gene_summary_path.stem,
                level=1,
                stage=Stage.ANALYSIS,
                method=Method.MAGECK,
                state=State.REPORT,
                ext="pdf",
            ).merge(tag)
        )

        output_dir = Path(outdir or self.default_output_dir(base_tag.root)).absolute()
        prefix = (
            output_prefix
            or f"{base_tag.root}.S{base_tag.level:02d}.{base_tag.stage}.{base_tag.method}.{base_tag.state}"
        )
        prefix_path = output_dir / prefix

        runner = self.run_plot(
            count_table=count_table_path,
            gene_summary=gene_summary_path,
            output_prefix=prefix_path,
            params=params,
            cfg=cfg,
        )

        key, name = self.build_element_name(base_tag, "plot")
        determinants = self.signature_determinants(params, subroutine="plot")

        return Element(
            key,
            runner,
            tag=base_tag,
            determinants=determinants,
            inputs=(count_table_path, gene_summary_path),
            artifacts={
                "output_prefix": prefix_path,
                "pdf": prefix_path.with_suffix(".pdf"),
            },
            pres=self._resolve_pres(count_table) + self._resolve_pres(gene_summary),
            name=name,
        )

    @subroutine
    def run_plot(
        self,
        count_table: Path | str,
        gene_summary: Path | str,
        output_prefix: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for ``mageck plot``."""
        count_table_path = Path(count_table).absolute()
        gene_summary_path = Path(gene_summary).absolute()
        prefix = Path(output_prefix).absolute()

        arguments = [
            "plot",
            "--count-table",
            str(count_table_path),
            "--gene-summary",
            str(gene_summary_path),
            "--output-prefix",
            str(prefix),
        ]
        arguments.extend(self.to_cli(params or Params(), subroutine="plot"))

        return (
            arguments,
            "plot",
            [count_table_path, gene_summary_path],
            [prefix.with_suffix(".pdf")],
            None,
            None,
            None,
        )

    @element
    def pathway(
        self,
        gene_ranking: Element | Path | str,
        gmt_file: Element | Path | str,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        output_prefix: str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> TableElement:
        """Run ``mageck pathway`` and return pathway summary table."""
        ranking_path = self._resolve_path(
            gene_ranking, preferred_keys=("gene_summary", "tsv", "path")
        )
        gmt_path = self._resolve_path(gmt_file, preferred_keys=("gmt", "path", "tsv"))

        base_tag = (
            from_prior(
                gene_ranking.tag,
                tag,
                stage=Stage.ANALYSIS,
                method=Method.MAGECK,
                state=State.ENRICHMENT,
                ext="txt",
            )
            if isinstance(gene_ranking, Element)
            else ElementTag(
                root=ranking_path.stem,
                level=1,
                stage=Stage.ANALYSIS,
                method=Method.MAGECK,
                state=State.ENRICHMENT,
                ext="txt",
            ).merge(tag)
        )

        output_dir = Path(outdir or self.default_output_dir(base_tag.root)).absolute()
        prefix = (
            output_prefix
            or f"{base_tag.root}.S{base_tag.level:02d}.{base_tag.stage}.{base_tag.method}.{base_tag.state}"
        )
        prefix_path = output_dir / prefix
        pathway_summary = prefix_path.with_suffix(".pathway_summary.txt")

        runner = self.run_pathway(
            gene_ranking=ranking_path,
            gmt_file=gmt_path,
            output_prefix=prefix_path,
            params=params,
            cfg=cfg,
        )

        key, name = self.build_element_name(base_tag, "pathway")
        determinants = self.signature_determinants(params, subroutine="pathway")

        return TableElement(
            key,
            runner,
            tag=base_tag,
            tsv=pathway_summary,
            artifacts={
                "tsv": pathway_summary,
                "pathway_summary": pathway_summary,
                "output_prefix": prefix_path,
            },
            determinants=determinants,
            inputs=(ranking_path, gmt_path),
            pres=self._resolve_pres(gene_ranking) + self._resolve_pres(gmt_file),
            name=name,
        )

    @subroutine
    def run_pathway(
        self,
        gene_ranking: Path | str,
        gmt_file: Path | str,
        output_prefix: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for ``mageck pathway``."""
        ranking_path = Path(gene_ranking).absolute()
        gmt_path = Path(gmt_file).absolute()
        prefix = Path(output_prefix).absolute()

        arguments = [
            "pathway",
            "--gene-ranking",
            str(ranking_path),
            "--gmt-file",
            str(gmt_path),
            "--output-prefix",
            str(prefix),
        ]
        arguments.extend(self.to_cli(params or Params(), subroutine="pathway"))

        return (
            arguments,
            "pathway",
            [ranking_path, gmt_path],
            [prefix.with_suffix(".pathway_summary.txt")],
            None,
            None,
            None,
        )


################################################################################
# Convenience functions for table juggling
################################################################################


def calculate_onesided_fdr(
    method: str = "fdr_bh", p_column: str = "p.high", column_name: str = "FDR High"
) -> Callable[[DataFrame], DataFrame]:
    """
    Calculate one-sided FDR from two-sided p-values using the Benjamini-Hochberg procedure.
    """

    def _calculate_fdr(df: DataFrame) -> DataFrame:
        df[column_name] = multipletests(df[p_column], method=method)[1]
        return df[[column_name]]

    return _calculate_fdr


def mageck_melt(
    metrics: list[str] = ["LFC", "Score", "p-value", "FDR"],
    set_name: str = "cell",
    additional_columns_per_set: (
        tuple[Iterable[str], Callable[[DataFrame], DataFrame]] | None
    ) = None,
) -> Callable[[DataFrame], DataFrame]:
    """
    Create a mageck input table from a long-format DataFrame.
    """

    def _melt(df: DataFrame) -> DataFrame:
        id_cols = ["sgRNA Number", "sgRNA", "Sensor", "Sequence"]

        value_cols = [
            c
            for c in df.columns
            if "(" in c and ")" in c and any(m in c for m in metrics)
        ]

        long_df = df.melt(
            id_vars=id_cols,
            value_vars=value_cols,
            var_name="variable",
            value_name="value",
        )
        long_df[["metric", set_name]] = long_df["variable"].str.extract(
            r"^(.*?) \((.*?)\)$"
        )
        long_df = long_df.drop(columns=["variable"])
        wide_df = long_df.pivot_table(
            index=id_cols + [set_name], columns="metric", values="value"
        ).reset_index()
        if additional_columns_per_set:
            wide_df[additional_columns_per_set[0]] = None
            for _, group_df in wide_df.groupby(set_name):
                additional_df = additional_columns_per_set[1](group_df)
                wide_df.loc[group_df.index, additional_columns_per_set[0]] = (
                    additional_df
                )
            assert (
                wide_df[additional_columns_per_set[0]].notna().all().all()
            ), "The additional column per set must be fully populated for all rows."

        return wide_df

    return _melt
