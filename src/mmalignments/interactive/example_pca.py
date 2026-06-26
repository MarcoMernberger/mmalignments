"""
example_pca.py
──────────────
Concrete usage example: interactive PCA plot exportable to standalone HTML.

Run interactively:
    python example_pca.py --serve

Export to HTML:
    python example_pca.py --export --outdir dist/
"""

from __future__ import annotations

import argparse

import numpy as np  # type: ignore[import]
import pandas as pd  # type: ignore[import]
import param  # type: ignore[import]
import plotly.express as px  # type: ignore[import]
from sklearn.decomposition import PCA  # type: ignore[import]
from sklearn.preprocessing import StandardScaler  # type: ignore[import]

from mmalignments.interactive.panels import (
    BaseAnalysis,
    BasePlot,
    PanelPlot,
    PlotState,
)

# ─────────────────────────────────────────────────────────────────────────────
# 1.  State  — all tunable parameters, grouped by prefix
# ─────────────────────────────────────────────────────────────────────────────


class PCAState(PlotState):
    # ── analysis params  (prefix: analysis_) ─────────────────────────────────
    analysis_n_components = param.Integer(
        default=2, bounds=(2, 10), doc="Number of PCA components to compute"
    )
    analysis_scale = param.Boolean(
        default=True, doc="Standardise features before PCA"
    )  # noqa: E501
    analysis_min_variance = param.Number(
        default=0.0,
        bounds=(0.0, 1.0),
        step=0.01,
        doc="Drop features with variance below this threshold",
    )

    # ── plot layout params  (prefix: plot_) ───────────────────────────────────
    plot_pc_x = param.Integer(
        default=1, bounds=(1, 10), doc="PC to show on X axis"
    )  # noqa: E501
    plot_pc_y = param.Integer(
        default=2, bounds=(1, 10), doc="PC to show on Y axis"
    )  # noqa: E501
    plot_color_by = param.Selector(
        default="genotype",
        objects=["genotype", "tissue", "replicate"],
        doc="Column used for point color",
    )
    plot_symbol_by = param.Selector(
        default="tissue",
        objects=["genotype", "tissue", "replicate", "none"],
        doc="Column used for point symbol",
    )
    plot_marker_size = param.Integer(default=10, bounds=(3, 30))
    plot_show_loadings = param.Boolean(
        default=False, doc="Overlay top-5 feature loadings as arrows"
    )
    plot_title = param.String(default="PCA")
    io_output_dir = param.String(
        default="outputs",
        doc="Directory for saving outputs (e.g. exported HTML)",
    )
    io_stem = param.String(
        default="pca",
        doc="Filename stem for outputs (e.g. exported HTML)",
    )


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Analysis  — pure computation, receives state as read-only namespace
# ─────────────────────────────────────────────────────────────────────────────


class PCAAnalysis(BaseAnalysis):
    """
    Runs PCA on a counts matrix (samples × features).
    Expects input df to have metadata columns + numeric feature columns.
    Returns a tidy df with PC1…PCn, explained variance, and metadata.
    """

    _META_COLS = ["sample", "genotype", "tissue", "replicate"]

    def run(self, data: pd.DataFrame, state: PCAState) -> pd.DataFrame:  # type: ignore[override]
        meta = data[self._META_COLS].copy()
        counts = data.drop(columns=self._META_COLS).astype(float)

        # variance filter
        if state.analysis_min_variance > 0:
            counts = counts.loc[:, counts.var() >= state.analysis_min_variance]

        X = (
            StandardScaler().fit_transform(counts)
            if state.analysis_scale
            else counts.values
        )

        n = min(state.analysis_n_components, X.shape[0], X.shape[1])
        pca = PCA(n_components=n)
        coords = pca.fit_transform(X)

        pc_cols = {f"PC{i+1}": coords[:, i] for i in range(n)}
        ev_pct = {
            f"PC{i+1}_var": round(pca.explained_variance_ratio_[i] * 100, 1)
            for i in range(n)
        }

        result = pd.concat(
            [meta, pd.DataFrame(pc_cols, index=data.index)], axis=1
        )  # noqa: E501
        result.attrs["explained_variance"] = ev_pct
        result.attrs["loadings"] = pd.DataFrame(
            pca.components_.T,
            index=counts.columns,
            columns=[f"PC{i+1}" for i in range(n)],
        )
        return result


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Plot  — pure rendering, receives tidy df + state
# ─────────────────────────────────────────────────────────────────────────────


class PCAPlot(BasePlot):

    def render(self, df: pd.DataFrame, state: PCAState):  # type: ignore[override]
        xc = f"PC{state.plot_pc_x}"
        yc = f"PC{state.plot_pc_y}"

        # graceful fallback if requested PC wasn't computed
        available = [
            c for c in df.columns if c.startswith("PC") and "_" not in c
        ]  # noqa: E501
        if xc not in available:
            xc = available[0]
        if yc not in available:
            yc = available[min(1, len(available) - 1)]

        ev = df.attrs.get("explained_variance", {})
        xlabel = f"{xc}  ({ev.get(xc + '_var', '?')}% var)"
        ylabel = f"{yc}  ({ev.get(yc + '_var', '?')}% var)"

        sym_col = (
            None if state.plot_symbol_by == "none" else state.plot_symbol_by
        )  # noqa: E501

        fig = px.scatter(
            df,
            x=xc,
            y=yc,
            color=state.plot_color_by,
            symbol=sym_col,
            hover_data=["sample", "genotype", "tissue"],
            title=state.plot_title,
            labels={xc: xlabel, yc: ylabel},
        )
        fig.update_traces(marker={"size": state.plot_marker_size})
        fig.update_layout(margin={"l": 40, "r": 20, "t": 50, "b": 40})

        # optional loading arrows
        if state.plot_show_loadings:
            loadings = df.attrs.get("loadings")
            if loadings is not None and xc in loadings and yc in loadings:
                top = loadings[[xc, yc]].abs().sum(axis=1).nlargest(5).index
                scale = df[[xc, yc]].abs().max().max() * 0.4
                for feat in top:
                    lx = float(loadings.loc[feat, xc]) * scale
                    ly = float(loadings.loc[feat, yc]) * scale
                    fig.add_annotation(
                        x=lx,
                        y=ly,
                        ax=0,
                        ay=0,
                        xref="x",
                        yref="y",
                        axref="x",
                        ayref="y",
                        text=feat,
                        showarrow=True,
                        arrowhead=2,
                        arrowcolor="gray",
                        font={"size": 9, "color": "gray"},
                    )
        return fig


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Demo data  (mimics your {genotype}_{tissue}_{replicate} format)
# ─────────────────────────────────────────────────────────────────────────────


def make_demo_data(n_features: int = 500, seed: int = 42) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    samples = (
        [
            {
                "sample": f"p53-R172H_NSCLC_{r}",
                "genotype": "p53-R172H",
                "tissue": "NSCLC",
                "replicate": str(r),
            }
            for r in range(1, 4)
        ]
        + [
            {
                "sample": f"WT_NSCLC_{r}",
                "genotype": "WT",
                "tissue": "NSCLC",
                "replicate": str(r),
            }
            for r in range(1, 4)
        ]
        + [
            {
                "sample": f"p53-R172H_Liver_{r}",
                "genotype": "p53-R172H",
                "tissue": "Liver",
                "replicate": str(r),
            }
            for r in range(1, 4)
        ]
        + [
            {
                "sample": f"WT_Liver_{r}",
                "genotype": "WT",
                "tissue": "Liver",
                "replicate": str(r),
            }
            for r in range(1, 4)
        ]
    )
    meta = pd.DataFrame(samples)
    n = len(meta)

    # add structured signal: genotype effect on PC1, tissue effect on PC2
    geno_signal = (meta["genotype"] == "p53-R172H").astype(float).values
    tissue_signal = (meta["tissue"] == "NSCLC").astype(float).values

    features = rng.standard_normal((n, n_features))
    features[:, :50] += (
        geno_signal[:, None] * 3
    )  # genotype separates on first 50 genes    # noqa: E501
    features[:, 50:100] += tissue_signal[:, None] * 3  # tissue on next 50

    feat_df = pd.DataFrame(
        features, columns=[f"Gene{i:04d}" for i in range(n_features)]
    )
    return pd.concat([meta, feat_df], axis=1)


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Entry point
# ─────────────────────────────────────────────────────────────────────────────


def build_app(data: pd.DataFrame, **kwargs) -> PanelPlot:
    state = PCAState()
    for key, value in kwargs.items():
        if hasattr(state, key):
            setattr(state, key, value)
    print(state)
    # if outdir:
    #     state.output_dir = outdir
    return PanelPlot(
        data=data,
        analysis=PCAAnalysis(),
        plot=PCAPlot(),
        state=state,
    )


# .from_state(
#         state_cls=PCAState,
#         data=data,
#     )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--serve", action="store_true", help="Launch Panel server")
    parser.add_argument(
        "--export", action="store_true", help="Export to standalone HTML"
    )
    parser.add_argument("--outdir", default="dist", help="Output directory for HTML")
    args = parser.parse_args()

    data = make_demo_data()

    if args.serve:
        app = build_app(data)
        app.show()
    elif args.export:
        app = build_app(data, io_output_dir=args.outdir, io_stem="pca")
        html = app.export_html()  # , , open_browser=True)
        print(f"\n✓ Standalone HTML: {html}")
    else:
        print("Use --serve or --export")
        print("Example: python example_pca.py --export --outdir dist/")
