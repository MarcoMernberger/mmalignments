from __future__ import annotations

from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt  # type: ignore[import]
import numpy as np
import pandas as pd  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

from mmalignments.models.elements import (
    Element,
    element,
    generate_element_key_name,
)
from mmalignments import __version__
from mmalignments.models.externals import Runnable
from mmalignments.models.parameters import Params
from mmalignments.models.tags import (
    ElementTag,
    Method,
    PartialElementTag,
    Stage,
    State,
    from_prior,
)
from mmalignments.services.io import from_json, write_json
from mmalignments.services.plots import (  # type: ignore[import]
    FIGURE_SUFFIXES,
    save_figure,
)


def _mouse_id_from_sample(sample: str) -> str:
    # Example: PRLC_126624_SCC_1 -> PRLC_126624
    parts = sample.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else sample


def _tissue_from_sample(sample: str) -> str:
    # Example: PRLC_126624_Kidney / PRLC_126624_SCC_1
    if "Kidney" in sample:
        return "Kidney"
    if "SCC" in sample:
        return "SCC"
    return "Other"


def _safe_float(x: Any) -> float | None:
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def _parse_picard_markdup_metrics(metrics_txt: str | Path) -> float | None:
    """
    Picard/GATK MarkDuplicates metrics: the duplication rate is usually in
    column PERCENT_DUPLICATION.
    We'll parse the first non-comment data line under the header.
    """
    p = Path(metrics_txt)
    if not p.exists():
        return None

    header = None
    for line in p.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        if (
            header is None
            and line.startswith("LIBRARY")
            and "PERCENT_DUPLICATION" in line
        ):
            header = line.split("\t")
            continue
        if header is not None and not line.startswith("##"):
            fields = line.split("\t")
            if len(fields) != len(header):
                continue
            d = dict(zip(header, fields))
            return _safe_float(d.get("PERCENT_DUPLICATION"))
    return None


def _parse_picard_hs_metrics(
    metrics_txt: str | Path,
) -> tuple[float | None, float | None]:
    """
    CollectHsMetrics has columns like PCT_SELECTED_BASES and
    MEAN_TARGET_COVERAGE. We'll parse similarly.
    Returns: (on_target_pct, mean_target_depth)
    """
    p = Path(metrics_txt)
    if not p.exists():
        return (None, None)

    header = None
    for line in p.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        # HS metrics header line often starts with BAIT_SET / or has PCT_SELECTED_BASES  # noqa: E501
        if header is None and (
            "PCT_SELECTED_BASES" in line and "MEAN_TARGET_COVERAGE" in line
        ):
            header = line.split("\t")
            continue
        if header is not None:
            fields = line.split("\t")
            if len(fields) != len(header):
                continue
            d = dict(zip(header, fields))
            on_target = _safe_float(d.get("PCT_SELECTED_BASES"))
            mean_depth = _safe_float(d.get("MEAN_TARGET_COVERAGE"))
            return (on_target, mean_depth)

    return (None, None)


@dataclass
class MutationalLoadReport:

    output_dir: Path = Path("results/report")
    name: str = "MutationalLoadReport"
    version: str = __version__

    ############################################################################
    # Helpers
    ############################################################################

    def build_element_name(
        self, tag: ElementTag, subcommand: str, suffix: str = "", **param_str
    ) -> tuple[str, str]:
        return generate_element_key_name(tag, self.name, self.version)

    def __get_counts(
        self, counts: dict
    ) -> tuple[float | None, float | None, int, int, int]:
        # core fields (from count json)

        # core fields (from count json)
        callable_mb = _safe_float(counts.get("callable_mb"))
        n_total = int(counts.get("total", 0))
        n_snv = int(counts.get("snps", 0))
        n_indel = int(counts.get("indels", 0))

        mut_load = _safe_float(counts.get("mutational_load"))
        # fallback if mutational_load not stored
        if callable_mb and callable_mb > 0:
            mut_load = n_total / callable_mb
        return mut_load, callable_mb, n_total, n_snv, n_indel

    def __get_duprate(
        self, markdup_by_sample: Mapping[str, Element] | None, tumor: str = ""
    ) -> float | None:
        dup_rate = None
        if markdup_by_sample:
            qc = markdup_by_sample.get(tumor, None)
            dup_rate = (
                _parse_picard_markdup_metrics(qc.metrics) if qc else None  # noqa: E501
            )
        return dup_rate

    def __get_hsmetrics(
        self, hs_by_sample: Mapping[str, Element] | None = None, tumor: str = ""
    ) -> tuple[float | None, float | None]:
        on_target_pct = None
        mean_target_depth = None
        if hs_by_sample:
            qc = hs_by_sample.get(tumor, None)
            if qc:
                on_target_pct, mean_target_depth = _parse_picard_hs_metrics(
                    qc.metrics
                )  #   noqa: E501
        return on_target_pct, mean_target_depth

    def __collect_prerequisites(
        self, counted_by_tumor, markdup_by_sample, hs_by_sample
    ) -> tuple[Element]:
        inputs = []
        for el in counted_by_tumor.values():
            inputs.append(el)
        if markdup_by_sample:
            for el in markdup_by_sample.values():
                inputs.append(el)
        if hs_by_sample:
            for el in hs_by_sample.values():
                inputs.append(el)
        return tuple(inputs)

    def __collect_inputs(
        self,
        prerequisites: tuple[Element],
    ) -> tuple[Path]:
        inputs = []
        for el in prerequisites:
            if el.output_files:
                for path in el.output_files:
                    inputs.append(path)
        return tuple(inputs)

    ############################################################################
    # Write tables
    ############################################################################

    def build_report_table(
        self,
        pairs: Mapping[str, str],  # tumor normal pairs
        counted_by_tumor: Mapping[str, Element],  # the final mutational loads
        *,
        markdup_by_sample: Mapping[str, Element] | None = None,
        hs_by_sample: Mapping[str, Element] | None = None,
    ) -> DataFrame:
        rows: list[dict] = []
        for tumor, normal in pairs.items():
            if tumor not in counted_by_tumor:
                raise ValueError(
                    f"Tumor sample {tumor} not found in counted_by_tumor"
                )  # noqa: E501

            count_elem = counted_by_tumor[tumor]
            counts_path = Path(count_elem.artifacts["json"])
            counts = from_json(counts_path)

            # optional QC
            mut_load, callable_mb, n_total, n_snv, n_indel = self.__get_counts(
                counts
            )  # noqa: E501
            dup_rate = self.__get_duprate(markdup_by_sample, tumor)
            on_target_pct, mean_target_depth = self.__get_hsmetrics(
                hs_by_sample, tumor
            )  # noqa: E501
            row = {
                "sample": tumor,
                "normal": normal,
                "mouse": _mouse_id_from_sample(tumor),
                "tissue": _tissue_from_sample(tumor),
                "callable_mb": callable_mb,
                "n_snv": n_snv,
                "n_indel": n_indel,
                "n_total": n_total,
                "mut_load": mut_load,
                "dup_rate": dup_rate,
                "on_target_pct": on_target_pct,
                "mean_target_depth": mean_target_depth,
            }
            rows.append(row)

        columns = [
            "sample",
            "normal",
            "mouse",
            "tissue",
            "callable_mb",
            "n_snv",
            "n_indel",
            "n_total",
            "mut_load",
            "dup_rate",
            "on_target_pct",
            "mean_target_depth",
        ]
        return DataFrame(rows, columns=columns)

    def write_tables(
        self,
        pairs: Mapping[str, str],  # tumor normal pairs
        counted_by_tumor: Mapping[str, Element],  # the final mutational loads
        out_tsv: Path,
        out_json: Path,
        *,
        markdup_by_sample: Mapping[str, Element] | None = None,
        hs_by_sample: Mapping[str, Element] | None = None,
    ) -> Runnable:
        def run():
            df = self.build_report_table(
                pairs,
                counted_by_tumor,
                markdup_by_sample=markdup_by_sample,
                hs_by_sample=hs_by_sample,
            )
            df.to_csv(out_tsv, sep="\t", index=False)
            dict_json = {
                "rows": df.to_dict(orient="records"),
                "columns": df.columns.tolist(),
            }
            write_json(out_json, dict_json)

        return Runnable(run, [], "write mutational load report")

    @element
    def write(
        self,
        pairs: Mapping[str, str],  # tumor normal pairs
        counted_by_tumor: Mapping[str, Element],  # the final mutational loads
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        markdup_by_sample: Mapping[str, Element] | None = None,
        hs_by_sample: Mapping[str, Element] | None = None,
        params: Params | None = None,
    ) -> Element:

        out_dir = Path(outdir or self.output_dir).absolute()
        tag = from_prior(
            next(iter(counted_by_tumor.values())).tag,
            tag,
            stage=Stage.SUMMARY,
            method=Method.CUSTOM,
            state=State.REPORT,
            ext="json",
        )
        outname = filename or tag.default_output
        out_json = (out_dir / outname).with_suffix(".json")
        out_tsv = out_json.with_suffix(".tsv")

        runner = self.write_tables(
            pairs,
            counted_by_tumor,
            out_tsv,
            out_json,
            markdup_by_sample=markdup_by_sample,
            hs_by_sample=hs_by_sample,
        )
        key, name = self.build_element_name(tag, "write")

        pres = self.__collect_prerequisites(
            counted_by_tumor, markdup_by_sample, hs_by_sample
        )
        inputs = self.__collect_inputs(pres)
        determinants = ()

        return Element(
            key=key,
            run=runner,
            tag=tag,
            determinants=determinants,
            inputs=inputs,
            artifacts={
                "tsv": out_tsv,
                "json": out_json,
            },
            pres=tuple(pres),
            name=name,
        )

    ############################################################################
    # Write Plots
    ############################################################################

    def plot_grouped(
        self,
        report_table: Path | str,
        out_plot: Path | str,
        *,
        group_col: str,
        value_col: str = "mut_load",
        mode: str = "bar",  # "bar" | "box"
        error: str = "std",  # "std" | "sem" | None
        params=None,
    ):
        def __run():
            suffixes = params.get("suffixes", None) if params else None

            df = pd.read_csv(report_table, sep="\t")

            # drop NaNs
            df = df[[group_col, value_col]].dropna()

            groups = df.groupby(group_col)[value_col].apply(list)

            labels = sorted(groups.index)
            data = [groups[ll] for ll in labels]

            f = plt.figure()

            if mode == "box":
                plt.boxplot(data, labels=labels, showfliers=True)

            elif mode == "bar":
                means = []
                errs = []

                for vals in data:
                    vals = np.array(vals, dtype=float)
                    means.append(vals.mean())

                    if len(vals) > 1 and error is not None:
                        if error == "std":
                            errs.append(vals.std())
                        elif error == "sem":
                            errs.append(vals.std() / np.sqrt(len(vals)))
                        else:
                            errs.append(0)
                    else:
                        errs.append(0)

                x = np.arange(len(labels))

                plt.bar(x, means, yerr=errs if any(errs) else None, capsize=5)
                plt.xticks(x, labels, rotation=45, ha="right")

            else:
                raise ValueError(f"Unknown mode: {mode}")

            plt.ylabel(f"{value_col} (variants / callable Mb)")
            plt.title(f"{value_col} grouped by {group_col}")
            plt.tight_layout()

            save_figure(f, Path(out_plot), suffixes=suffixes)
            plt.close()

        return Runnable(__run, [], f"plot {value_col} by {group_col}")

    # def plot_by_group(
    #     self,
    #     report_table: Path | str,
    #     out_plot_group: Path | str,
    #     *,
    #     params: Params | None = None,
    # ) -> Runnable:
    #     # ---------- plots ----------
    #     # Plot 1: by group (Kidney vs SCC vs Other)
    #     # (Even if you only have SCC tumors, it will still plot.)
    #     def __run():
    #         suffixes = params.get("suffixes", None) if params else None
    #         groups = {}
    #         df = pd.read_csv(report_table, sep="\t")
    #         for _, r in df.iterrows():
    #             g = r["tissue"]
    #             v = r["mut_load"]
    #             if v is None:
    #                 continue
    #             groups.setdefault(g, []).append(float(v))

    #         if groups:
    #             labels = sorted(
    #                 groups.keys(),
    #                 key=lambda x: (
    #                     ["Kidney", "SCC", "Other"].index(x)
    #                     if x in ["Kidney", "SCC", "Other"]
    #                     else 99
    #                 ),
    #             )
    #             data = [groups[ll] for ll in labels]

    #             f = plt.figure()
    #             plt.boxplot(data, labels=labels, showfliers=True)
    #             plt.ylabel("Mutational load (variants / callable Mb)")
    #             plt.title("Mutational load by tissue group")
    #             plt.tight_layout()
    #             plt.close()
    #             save_figure(f, Path(out_plot_group), suffixes=suffixes)

    #     return Runnable(__run, [], "plot mutational load by group")

    # def plot_by_mouse(
    #     self,
    #     report_table: Path | str,
    #     out_plot_mouse: Path | str,
    #     *,
    #     params: Params | None = None,
    # ) -> Runnable:
    #     # Plot 2: by mouse (median per mouse as bar plot)
    #     print(out_plot_mouse)
    #     # raise ValueError("Intentional error to test dry run")  # noqa: E501

    #     def __run():
    #         suffixes = params.get("suffixes", None) if params else None
    #         df = pd.read_csv(report_table, sep="\t")
    #         mouse_vals = {}
    #         for _, r in df.iterrows():
    #             mid = r["mouse_id"]
    #             v = r["mut_load"]
    #             if v is None:
    #                 continue
    #             mouse_vals.setdefault(mid, []).append(float(v))

    #         if mouse_vals:
    #             mouse_ids = sorted(mouse_vals.keys())
    #             # median
    #             meds = []
    #             for mid in mouse_ids:
    #                 xs = sorted(mouse_vals[mid])
    #                 n = len(xs)
    #                 meds.append(
    #                     xs[n // 2]
    #                     if n % 2 == 1
    #                     else (xs[n // 2 - 1] + xs[n // 2]) / 2.0
    #                 )

    #             f = plt.figure()
    #             plt.bar(mouse_ids, meds)
    #             plt.xticks(rotation=45, ha="right")
    #             plt.ylabel("Median mutational load (variants / callable Mb)")
    #             plt.title("Mutational load by mouse (median across tumors)")
    #             plt.tight_layout()
    #             print(Path(out_plot_mouse))
    #             save_figure(f, Path(out_plot_mouse), suffixes=suffixes)
    #             plt.close()

    #     return Runnable(__run, [], "plot mutational load by mouse")

    @element
    def bygroup(
        self,
        report_table: Element,
        *,
        group: str = "tissue",
        value: str = "mut_load",
        mode: str = "bar",
        error: str = "std",
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
    ) -> Element:
        out_dir = Path(outdir or self.output_dir).absolute()
        param = tag.param if tag and hasattr(tag, "param") else f"by_{group}"
        tag = from_prior(
            report_table.tag,
            tag,
            stage=Stage.SUMMARY,
            method=Method.CUSTOM,
            state=State.REPORT,
            param=param,
            ext="png",
        )
        out_plot_group = out_dir / tag.default_output
        if filename:
            filename = Path(filename)
            out_plot_group = filename.with_suffix(f".{tag.param}{filename.suffix}")

        runner = self.plot_grouped(
            report_table.artifacts["tsv"],
            out_plot_group,
            group_col=group,
            value_col=value,
            mode=mode,
            error=error,
            params=params,
        )

        key, name = self.build_element_name(tag, "bygroup")

        inputs = (report_table.artifacts["tsv"],)
        determinants = ()
        suffixes = (
            params.get("suffixes", FIGURE_SUFFIXES)
            if params
            else FIGURE_SUFFIXES  # noqa: E501
        )
        artifacts = {}
        for suffix in suffixes:
            artifacts[suffix.lstrip(".")] = out_plot_group.with_suffix(suffix)
        return Element(
            key=key,
            run=runner,
            tag=tag,
            determinants=determinants,
            inputs=inputs,
            artifacts=artifacts,
            pres=(report_table,),
            name=name,
        )

    # @element
    # def bymouse(
    #     self,
    #     report_table: Element,
    #     *,
    #     tag: PartialElementTag | ElementTag | None = None,
    #     outdir: Path | str | None = None,
    #     filename: Path | str | None = None,
    #     params: Params | None = None,
    # ) -> Element:
    #     out_dir = Path(outdir or self.output_dir).absolute()
    #     tag = from_prior(
    #         report_table.tag,
    #         tag,
    #         stage=Stage.SUMMARY,
    #         method=Method.CUSTOM,
    #         state=State.REPORT,
    #         param="bymouse",
    #         ext="png",
    #     )

    #     out_plot_mouse = out_dir / tag.default_output
    #     if filename:
    #         filename = Path(filename)
    #         out_plot_mouse = filename.with_suffix(f".{tag.param}{filename.suffix}")

    #     runner = self.plot_by_mouse(
    #         report_table.artifacts["tsv"], out_plot_mouse, params=params
    #     )
    #     key, name = self.build_element_name(tag, "bymouse")

    #     inputs = (report_table.artifacts["tsv"],)
    #     determinants = ()
    #     suffixes = (
    #         params.get("suffixes", FIGURE_SUFFIXES)
    #         if params
    #         else FIGURE_SUFFIXES  # noqa: E501
    #     )
    #     artifacts = {}
    #     for suffix in suffixes:
    #         artifacts[suffix.lstrip(".")] = out_plot_mouse.with_suffix(suffix)
    #     return Element(
    #         key=key,
    #         run=runner,
    #         tag=tag,
    #         determinants=determinants,
    #         inputs=inputs,
    #         artifacts=artifacts,
    #         pres=(report_table,),
    #         name=name,
    #     )

    # def plot_by_mouse(
    #     self,
    #     report_table: Path | str,
    #     out_plot_mouse: Path | str,
    #     *,
    #     params: Params | None = None,
    # ) -> Runnable:
    #     # Plot 2: by mouse (median per mouse as bar plot)
    #     print(out_plot_mouse)
    #     # raise ValueError("Intentional error to test dry run")  # noqa: E501

    #     def __run():
    #         suffixes = params.get("suffixes", None) if params else None
    #         df = pd.read_csv(report_table, sep="\t")
    #         mouse_vals = {}
    #         for _, r in df.iterrows():
    #             mid = r["mouse_id"]
    #             v = r["mut_load"]
    #             if v is None:
    #                 continue
    #             mouse_vals.setdefault(mid, []).append(float(v))

    #         if mouse_vals:
    #             mouse_ids = sorted(mouse_vals.keys())
    #             # median
    #             meds = []
    #             for mid in mouse_ids:
    #                 xs = sorted(mouse_vals[mid])
    #                 n = len(xs)
    #                 meds.append(
    #                     xs[n // 2]
    #                     if n % 2 == 1
    #                     else (xs[n // 2 - 1] + xs[n // 2]) / 2.0
    #                 )

    #             f = plt.figure()
    #             plt.bar(mouse_ids, meds)
    #             plt.xticks(rotation=45, ha="right")
    #             plt.ylabel("Median mutational load (variants / callable Mb)")
    #             plt.title("Mutational load by mouse (median across tumors)")
    #             plt.tight_layout()
    #             print(Path(out_plot_mouse))
    #             save_figure(f, Path(out_plot_mouse), suffixes=suffixes)
    #             plt.close()

    #     return Runnable(__run, [], "plot mutational load by mouse")

    ############################################################################
    # Build report
    ############################################################################

    def build(
        self,
        pairs: Mapping[str, str],  # tumor normal pairs
        counted_by_tumor: Mapping[str, Element],  # the final mutational loads
        *,
        markdup_by_sample: Mapping[str, Element] | None = None,
        hs_by_sample: Mapping[str, Element] | None = None,
        tags: Mapping[str, PartialElementTag | ElementTag] | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        parameters: Mapping[str, Params] | None = None,
    ) -> tuple[Element, ...]:
        out_dir = outdir or Path(self.output_dir).absolute()

        # write the report table first, then use it for plotting
        report_write = self.write(
            pairs,
            counted_by_tumor,
            tag=tags.get("write") if tags else None,
            outdir=out_dir,
            filename=filename,
            markdup_by_sample=markdup_by_sample,
            hs_by_sample=hs_by_sample,
            params=parameters.get("write") if parameters else None,
        )
        # TODO: make generic
        # plot by tissue (default)
        report_plot_group = self.bygroup(
            report_write,
            outdir=out_dir,
            filename=filename,
            params=parameters.get("tissue") if parameters else None,
        )
        # plot by mouse
        report_plot_mouse = self.bygroup(
            report_write,
            group="mouse",
            outdir=out_dir,
            filename=filename,
            params=parameters.get("mouse") if parameters else None,
        )
        # plot by sample
        report_plot_sample = self.bygroup(
            report_write,
            group="sample",
            outdir=out_dir,
            filename=filename,
            params=parameters.get("sample") if parameters else None,
        )
        return report_write, report_plot_group, report_plot_mouse, report_plot_sample
