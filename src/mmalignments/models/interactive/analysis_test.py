from spec import PlotConfig, PlotLayer
from export import export_standalone_script
from spec import Analysis, AnalysisResult, PlotConfig, PlotLayer
from pathlib import Path
import pandas as pd
import numpy as np

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
    print("data", data, sep="\n")
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
        agg = enrichment.reset_index()
        agg.columns = ["sgRNA", "Log2 Relative Enrichment Score"]
        agg["treatment"] = "Rezatapopt vs DMSO"  # for the color role, if used

        return AnalysisResult(raw=raw, agg=agg)


if __name__ == "__main__":
    data = make_demo_data()
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
        labels={"title": "CRISPR Screen QC"},
    )
    # _app = build_app(_data, config)
    # _app.panel().servable(title="CRISPR Screen QC")
    project_dir = Path(
        "/clara/ffs/e/20260519_AG_Stiewe_Imke_Bullwinkel_Rezatapopt_screen"
    )
    extra_imports = (
        "import numpy as np",
        "from pathlib import Path",
        "import pandas as pd",
    )
    export_standalone_script(
        data=data,
        config=config,
        analysis_cls=CrisprScreenAnalysis,
        outdir=project_dir / "standalone",
        stem="standalone_app",
        extra_imports=extra_imports,
        extra_plot_classes=(),
        title="Interactive Standalone Analysis",
    )
