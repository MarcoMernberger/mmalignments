"""
demo_crispr.py
────────────────
Reproduces the sgRNA count/frequency/enrichment example from the draft,
fully wired and runnable. Demonstrates:
  - two layers of the SAME plot_type ("rate") with different aliases
    -> proves namespacing-by-alias avoids collisions
  - one layer of a different plot_type ("hbar") reading a DIFFERENT
    AnalysisResult source ("agg" instead of "raw")
"""

import numpy as np
import pandas as pd

from spec import Analysis, AnalysisResult, PlotConfig, PlotLayer
from app import build_app

# ─────────────────────────────────────────────────────────────────────────────
# Demo data
# ─────────────────────────────────────────────────────────────────────────────


def make_demo_data(n_sgrnas=20, seed=7) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    sgRNAs = [f"sg{i:03d}" for i in range(n_sgrnas)]
    rows = []
    for sg in sgRNAs:
        for treatment in ["Rezatapopt", "DMSO"]:
            count = rng.integers(50, 5000)
            rows.append({"sgRNA": sg, "treatment": treatment, "Count": count})
    data = pd.DataFrame(rows)
    return data


# ─────────────────────────────────────────────────────────────────────────────
# Analysis: computes Frequency (raw, per-row) + enrichment (agg, per-sgRNA)
# ─────────────────────────────────────────────────────────────────────────────


class CrisprScreenAnalysis(Analysis):
    def run(self, data: pd.DataFrame, state) -> AnalysisResult:
        raw = data.copy()
        total_per_treatment = raw.groupby("treatment")["Count"].transform("sum")
        raw["Frequency"] = raw["Count"] / total_per_treatment
        pivot = raw.pivot(index="sgRNA", columns="treatment", values="Frequency")
        enrichment = np.log2(pivot["Rezatapopt"] / pivot["DMSO"])
        # raw["Log2 Relative Enrichment Score"] = raw["sgRNA"].map(enrichment)
        agg = enrichment.reset_index()
        agg.columns = ["sgRNA", "Log2 Relative Enrichment Score"]
        agg["treatment"] = "Rezatapopt vs DMSO"  # for the color role, if used
        agg = agg.sort_values("Log2 Relative Enrichment Score", ascending=False)
        print(
            "CrisprScreenAnalysis.run(): raw.shape =",
            raw.shape,
            "agg.shape =",
            agg.shape,
        )
        print("CrisprScreenAnalysis.run(): raw.head():\n", raw.head(10))
        print("CrisprScreenAnalysis.run(): agg.head():\n", agg.head(10))
        return AnalysisResult(raw=raw, agg=agg)


# ─────────────────────────────────────────────────────────────────────────────
# Config: count + frequency (same plot_type "rate", different alias) + hbar
# ─────────────────────────────────────────────────────────────────────────────

config = PlotConfig(
    analysis=CrisprScreenAnalysis(),
    layers=[
        PlotLayer(
            plot_type="rate",
            alias="count",
            roles={"x": "sgRNA", "y": "Count", "color": "treatment"},
            source="raw",
        ),
        PlotLayer(
            plot_type="rate",
            alias="frequency",
            roles={"x": "sgRNA", "y": "Frequency", "color": "treatment"},
            source="raw",
        ),
        PlotLayer(
            plot_type="hbar",
            alias="enrichment",
            roles={"x": "Log2 Relative Enrichment Score", "y": "sgRNA"},
            source="agg",
        ),
    ],
    select_param="analysis_layer_select",  # optional, only needed for SelectiveRenderer
    labels={"title": "CRISPR Screen QC"},
)


if __name__ == "__main__":
    data = make_demo_data()
    app = build_app(data, config)

    print("Combined state params (proves no collision between count_/frequency_):")
    for p in app.state_cls.param:
        if p != "name":
            print(" ", p)

    pane = app._figure_pane  # property, triggers pn.bind -> _render once
    print("\n✓ Rendered pane object:", type(pane))

# ── Module-level servable: lets `panel serve demo_crispr.py` find it ─────────
# Runs unconditionally at import time (not inside __main__), because
# `panel serve` imports this file as a module rather than executing it
# as a script.
_data = make_demo_data()
_app = build_app(_data, config)
_app.panel().servable(title="CRISPR Screen QC")
