from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt

from mmalignments.models.tasks import Element


def _read_json(p: Union[str, Path]) -> dict:
    with open(Path(p), "r") as f:
        return json.load(f)


def _write_tsv(rows: List[dict], out_tsv: Path, columns: List[str]) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write(
                "\t".join("" if r.get(c) is None else str(r.get(c)) for c in columns)
                + "\n"
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


def _safe_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def _parse_picard_markdup_metrics(metrics_txt: Union[str, Path]) -> Optional[float]:
    """
    Picard/GATK MarkDuplicates metrics: the duplication rate is usually in column PERCENT_DUPLICATION.
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
    metrics_txt: Union[str, Path],
) -> Tuple[Optional[float], Optional[float]]:
    """
    CollectHsMetrics has columns like PCT_SELECTED_BASES and MEAN_TARGET_COVERAGE.
    We'll parse similarly.
    Returns: (on_target_pct, mean_target_depth)
    """
    p = Path(metrics_txt)
    if not p.exists():
        return (None, None)

    header = None
    for line in p.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        # HS metrics header line often starts with BAIT_SET / or has PCT_SELECTED_BASES
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
    output_dir: Union[str, Path] = "results/report"

    def build(
        self,
        pairs: Dict[str, str],
        counted_by_tumor: Dict[str, Element],
        qc_by_sample: Optional[Dict[str, Dict[str, Union[str, Path]]]] = None,
    ) -> Element:
        out_dir = Path(self.output_dir).absolute()
        out_tsv = out_dir / "mutational_load_summary.tsv"
        out_json = out_dir / "mutational_load_summary.json"

        out_plot_group_png = out_dir / "mut_load_by_group.png"
        out_plot_group_pdf = out_dir / "mut_load_by_group.pdf"
        out_plot_mouse_png = out_dir / "mut_load_by_mouse.png"
        out_plot_mouse_pdf = out_dir / "mut_load_by_mouse.pdf"

        qc_by_sample = qc_by_sample or {}

        # Build prerequisites: all count elements
        pres = list(counted_by_tumor.values())

        def _run() -> None:
            out_dir.mkdir(parents=True, exist_ok=True)

            rows: List[dict] = []
            for tumor, normal in pairs.items():
                if tumor not in counted_by_tumor:
                    # If you want, raise here; I keep it robust.
                    continue

                count_elem = counted_by_tumor[tumor]
                counts_path = Path(count_elem.artifacts["json"])
                counts = _read_json(counts_path)

                # core fields (from your json)
                callable_mb = _safe_float(counts.get("callable_mb"))
                n_total = int(counts.get("n_total", 0))
                n_snv = int(counts.get("n_snv", 0))
                n_indel = int(counts.get("n_indel", 0))
                mut_load = _safe_float(counts.get("mutational_load"))

                # fallback if mutational_load not stored
                if mut_load is None and callable_mb and callable_mb > 0:
                    mut_load = n_total / callable_mb

                # optional QC
                qc = qc_by_sample.get(tumor, {})
                dup_rate = (
                    _parse_picard_markdup_metrics(qc["markdup_metrics"])
                    if "markdup_metrics" in qc
                    else None
                )
                on_target_pct, mean_target_depth = (
                    _parse_picard_hs_metrics(qc["hs_metrics"])
                    if "hs_metrics" in qc
                    else (None, None)
                )

                row = {
                    "sample": tumor,
                    "normal": normal,
                    "mouse_id": _mouse_id_from_sample(tumor),
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
                "mouse_id",
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
            _write_tsv(rows, out_tsv, columns)

            with open(out_json, "w") as f:
                json.dump({"rows": rows, "columns": columns}, f, indent=2)

            # ---------- plots ----------
            # Plot 1: by group (Kidney vs SCC vs Other)
            # (Even if you only have SCC tumors, it will still plot.)
            groups = {}
            for r in rows:
                g = r["tissue"]
                v = r["mut_load"]
                if v is None:
                    continue
                groups.setdefault(g, []).append(float(v))

            if groups:
                labels = sorted(
                    groups.keys(),
                    key=lambda x: (
                        ["Kidney", "SCC", "Other"].index(x)
                        if x in ["Kidney", "SCC", "Other"]
                        else 99
                    ),
                )
                data = [groups[ll] for ll in labels]

                plt.figure()
                plt.boxplot(data, labels=labels, showfliers=True)
                plt.ylabel("Mutational load (variants / callable Mb)")
                plt.title("Mutational load by tissue group")
                plt.tight_layout()
                plt.savefig(out_plot_group_png, dpi=200)
                plt.savefig(out_plot_group_pdf)
                plt.close()

            # Plot 2: by mouse (median per mouse as bar plot)
            mouse_vals = {}
            for r in rows:
                mid = r["mouse_id"]
                v = r["mut_load"]
                if v is None:
                    continue
                mouse_vals.setdefault(mid, []).append(float(v))

            if mouse_vals:
                mouse_ids = sorted(mouse_vals.keys())
                # median
                meds = []
                for mid in mouse_ids:
                    xs = sorted(mouse_vals[mid])
                    n = len(xs)
                    meds.append(
                        xs[n // 2]
                        if n % 2 == 1
                        else (xs[n // 2 - 1] + xs[n // 2]) / 2.0
                    )

                plt.figure()
                plt.bar(mouse_ids, meds)
                plt.xticks(rotation=45, ha="right")
                plt.ylabel("Median mutational load (variants / callable Mb)")
                plt.title("Mutational load by mouse (median across tumors)")
                plt.tight_layout()
                plt.savefig(out_plot_mouse_png, dpi=200)
                plt.savefig(out_plot_mouse_pdf)
                plt.close()

        return Element(
            name="final_mutational_load_report",
            key="final_mutational_load_report",
            run=_run,
            parameters={},
            input_files=[],
            artifacts={
                "tsv": out_tsv,
                "json": out_json,
                "plot_group_png": out_plot_group_png,
                "plot_group_pdf": out_plot_group_pdf,
                "plot_mouse_png": out_plot_mouse_png,
                "plot_mouse_pdf": out_plot_mouse_pdf,
            },
            pres=pres,
            store_attributes={"multi_core": False},
        )
