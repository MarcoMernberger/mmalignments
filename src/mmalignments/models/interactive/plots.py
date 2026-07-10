"""
plots.py
─────────
Concrete plot types + their PlotState classes.

Pattern per plot type:
  XxxState(PlotState)  → layout knobs ONLY (marker_size, linewidth, colors...)
  XxxPlot(BasePlot)    → ROLE_SPECS (data contract) + render()

TraceProcessor is shared post-processing: anything that's just "set a
plotly property from a state value" lives here ONCE instead of being
duplicated in every plot class's render().
"""

from __future__ import annotations

import bisect
import heapq
from typing import Any

import numpy as np
import pandas as pd
import param  # type: ignore[import]
import plotly.colors
import plotly.express as px
import plotly.graph_objects as go
from pandas import DataFrame, Series
from param import Parameterized  # type: ignore[import]
from textwrap import shorten
from .spec import (
    BasePlot,
    DynamicPlotState,
    DynamicSelectorSpec,
    PlotLayer,
    PlotState,
    RoleSpec,
)

from .actions import action_export_figure, action_export_data, action_export_app

# ─────────────────────────────────────────────────────────────────────────────
# TraceProcessor — shared appearance logic, gated by trace_capabilities
# ─────────────────────────────────────────────────────────────────────────────


class TraceProcessor:
    """
    Applies generic appearance settings to a figure.

    apply()        : per-LAYER marker properties (size, line width/color) —
                      called once per layer, from that layer's own state.
    apply_shared()  : figure-level properties (title, axis labels, tick
                      relabeling) — called ONCE after all layers are merged,
                      reading from the combined state directly.
    """

    def apply(self, fig: go.Figure, state, capabilities: set[str]) -> go.Figure:
        return self._update_traces(fig, state, capabilities)

    def apply_shared(
        self, fig: go.Figure, combined_state, layers: list[PlotLayer], result
    ) -> go.Figure:
        fig = self._apply_labels(fig, combined_state)
        fig = self._update_tick_labels(fig, combined_state, layers, result)
        return fig

    def _update_traces(self, fig, state, capabilities) -> go.Figure:
        marker_kwargs = {}
        if "marker.size" in capabilities and hasattr(state, "marker_size"):
            marker_kwargs["size"] = state.marker_size
        if "marker.line.width" in capabilities and hasattr(state, "linewidth"):
            marker_kwargs.setdefault("line", {})["width"] = state.linewidth
        if "marker.line.color" in capabilities and hasattr(state, "edgecolor"):
            marker_kwargs.setdefault("line", {})["color"] = state.edgecolor
        if marker_kwargs:
            fig.update_traces(marker=marker_kwargs)
        return fig

    def _apply_labels(self, fig, combined_state) -> go.Figure:
        if getattr(combined_state, "Appearance_appearance_xlabel", None):
            fig.update_xaxes(title=combined_state.Appearance_appearance_xlabel)
        if getattr(combined_state, "Appearance_appearance_ylabel", None):
            fig.update_yaxes(title=combined_state.Appearance_appearance_ylabel)
        if getattr(combined_state, "Appearance_appearance_title", None):
            fig.update_layout(title=combined_state.Appearance_appearance_title)
        return fig

    def _update_tick_labels(
        self, fig, combined_state, layers: list[PlotLayer], result
    ) -> go.Figure:
        """Relabel x-ticks using a different column than the one plotted
        (e.g. plot numeric `position`, show `sample` as tick text). Uses
        the FIRST layer's source frame as the reference for tick values."""
        xtick_col = getattr(combined_state, "Appearance_appearance_xtick", None)
        if xtick_col is None or not layers:
            return fig
        df = result[layers[0].source]
        if xtick_col not in df.columns or "position" not in df.columns:
            return fig
        first_group = df.drop_duplicates(subset="position").sort_values("position")
        fig.update_layout(
            xaxis={
                "tickmode": "array",
                "tickvals": first_group["position"].astype(int).tolist(),
                "ticktext": first_group[xtick_col].astype(str).tolist(),
                "automargin": True,
            }
        )
        return fig


# ─────────────────────────────────────────────────────────────────────────────
# RatePlot (scatter-style) — layout state + plot class
# ─────────────────────────────────────────────────────────────────────────────


class RatePlotState(PlotState):
    """Layout knobs for RatePlot. NOT data roles — those are in PlotLayer."""

    sidebar_group = "Layer"  # overridden per-alias anyway via namespacing

    marker_size = param.Integer(default=8, bounds=(1, 50), label="Marker Size")
    linewidth = param.Number(default=1.0, bounds=(0.0, 10.0), label="Marker Edge Width")
    edgecolor = param.Color(default="#000000", label="Marker Edge Color")
    opacity = param.Number(default=1.0, bounds=(0.0, 1.0), label="Opacity")
    section = None


class RatePlot(BasePlot):
    ROLE_SPECS = {
        "x": RoleSpec(name="x", dtype="categorical", description="X-axis variable"),
        "y": RoleSpec(name="y", dtype="continuous", description="Y-axis variable"),
        "color": RoleSpec(
            name="color",
            dtype="categorical",
            required=False,
            description="Color grouping variable",
        ),
    }
    state_cls = RatePlotState
    trace_capabilities = {"marker.size", "marker.line.width", "marker.line.color"}

    def render(
        self,
        df: pd.DataFrame,
        layer: PlotLayer,
        state: RatePlotState,
        processor: TraceProcessor,
        fig=None,
    ) -> go.Figure:
        roles = layer.roles
        kwargs = {}
        if "color" in roles:
            kwargs["color"] = roles["color"]
            color_map = layer.color_maps.get("color")
            if color_map:
                kwargs["color_discrete_map"] = color_map
        if layer.hover:
            kwargs["hover_data"] = layer.hover

        new_fig = px.scatter(df, x=roles["x"], y=roles["y"], **kwargs)
        new_fig.update_traces(opacity=state.opacity)
        new_fig = processor.apply(new_fig, state, self.trace_capabilities)
        return _merge(fig, new_fig)


# ─────────────────────────────────────────────────────────────────────────────
# GenePlot (very situational) — layout state + plot class
# ─────────────────────────────────────────────────────────────────────────────


class GenePlotState(PlotState):
    """Layout knobs for GenePlot. NOT data roles — those are in PlotLayer."""

    sidebar_group = "Layer"  # overridden per-alias anyway via namespacing
    track_coordinates = param.String(
        default="",
        label="Coordinates",
        doc="Genomic coordinates to display. this can be either a range separated by a '-' (e.g. `1000-2000`) or name of a feature (e.g. `Exon 1`, `Intron 3`, `TP53` or `sgRNA N76`) or a Type (e.g. sgRNA).",
    )
    track_spacing = param.Number(
        default=0.4,
        bounds=(0.0, 2.0),
        label="Spacing fo sgRNAs",
        doc="Vertical spacing between sgRNA tracks. This is a relative value, where 1.0 corresponds to the height of a gene track.",
    )

    track_head_x = param.Number(
        default=0.4,
        bounds=(0.0, 1.0),
        label="Relative arrow Head Width",
        doc="Relative width of the arrow head for sgRNA tracks.",
    )
    track_head_y = param.Number(
        default=0.2,
        bounds=(0.0, 1.0),
        label="Relative arrow Head Height",
        doc="Relative height of the arrow head for sgRNA tracks.",
    )
    # tracks_arrow_width = param.Number(
    #     default=1.0, bounds=(0.0, 10.0), label="Arrow Shaft Width"
    # )
    # edgecolor = param.Color(default="#000000", label="Marker Edge Color")
    track_opacity = param.Number(default=1.0, bounds=(0.0, 1.0), label="Opacity")
    track_metric = param.Selector(
        default="Score",
        objects=["Score", "LFC", "FDR High", "FDR Two-sided"],
        label="Metric",
        doc="Metric to use for sgRNA coloring.",
    )
    track_metric_colorscale = param.Selector(
        default="Plasma",
        objects=["RdYlGn", "Viridis", "RdBu", "Plasma", "Turbo", "Cividis"],
        label="sgRNA Colorscale",
        doc="Colorscale to use for sgRNA coloring.",
    )
    section = None


class GenePlot(BasePlot):
    ROLE_SPECS = {
        "start": RoleSpec(
            name="start",
            dtype="continuous",
            required=False,
            description="Genomic start position column",
        ),
        "stop": RoleSpec(
            name="stop",
            dtype="continuous",
            required=False,
            description="Genomic stop position column",
        ),
        "seq": RoleSpec(
            name="seq",
            dtype="categorical",
            required=False,
            description="Sequence column (for hover and text labels)",
        ),
        "track": RoleSpec(
            name="track", dtype="categorical", description="Track column"
        ),
        "feature": RoleSpec(
            name="feature",
            dtype="categorical",
            required=False,
            description="Feature column (for hover and text labels)",
        ),
        "dir": RoleSpec(
            name="dir",
            dtype="categorical",
            required=False,
            description="Direction column (for hover and text labels)",
        ),
        "color": RoleSpec(
            name="color",
            dtype="continuous",
            required=False,
            description="Metric used to color sgRNA markers",
        ),
        "track_offset": RoleSpec(
            name="track_offset",
            dtype="continuous",
            required=False,
            description="Offset for sgRNA markers to avoid overlap",
        ),
        "ensembl": RoleSpec(
            name="ensembl",
            dtype="categorical",
            required=False,
            description="Ensembl ID column (for hover and text labels)",
        ),
    }
    state_cls = GenePlotState
    # trace_capabilities = {"marker.size", "marker.line.width", "marker.line.color"}

    def render(
        self,
        df: pd.DataFrame,
        layer: PlotLayer,
        state: GenePlotState,
        processor: TraceProcessor,
        fig=None,
    ) -> go.Figure:
        roles = layer.roles

        if (
            roles["start"] not in df.columns
            or roles["stop"] not in df.columns
            or roles["track"] not in df.columns
        ):
            return _merge(fig, go.Figure())

        y_role_col = roles.get("y")
        color_col = roles.get("color") if "color" in roles else None

        work = df.copy()
        work[roles["track"]] = work[roles["track"]].astype(str)  # .str.lower()
        xmin = (
            work.loc[work[roles["start"]] > 0, [roles["start"], roles["stop"]]]
            .min()
            .min()
        )
        xmax = (
            work.loc[work[roles["start"]] > 0, [roles["start"], roles["stop"]]]
            .max()
            .max()
        )
        # Map raw feature types to browser-like tracks.
        track_map = {
            "gene": "Gene",
            "transcript": "Transcript",
            "intron": "Transcript",
            "exon": "Transcript",
            "sequence": "Sequence",
            "sgRNA": "sgRNA",
        }
        if y_role_col and y_role_col in work.columns:
            work["_track"] = work[y_role_col].astype(str)
        else:
            work["_track"] = work[roles["track"]].map(track_map).fillna("Other")
        work["_track"] = work.apply(
            lambda r: track_map.get(r[roles["track"]], r["_track"]), axis=1
        )
        preferred = ["Gene", "Transcript", "Sequence", "sgRNA"]
        discovered = [
            t for t in work["_track"].dropna().unique().tolist() if t not in preferred
        ]
        track_order = [t for t in preferred if t in work["_track"].values] + discovered
        if not track_order:
            track_order = preferred

        n_tracks = len(track_order) + np.ceil(
            work[roles["track_offset"]].max() * state.track_spacing
        )  # add extra space for any unknown tracks
        track_y = {name: float(n_tracks - idx) for idx, name in enumerate(track_order)}
        new_fig = go.Figure()

        # Gene and transcript structures are rendered as interval shapes.
        for _, row in work.iterrows():
            ftype = row[roles["track"]]
            try:
                start = int(row[roles["start"]])
                stop = int(row[roles["stop"]])
            except Exception:
                continue

            if ftype == "gene":
                y = track_y.get("Gene")
                if y is None:
                    continue

                new_fig.add_trace(
                    go.Scatter(
                        x=[
                            start,
                            start,
                            stop,
                            stop,
                        ],
                        y=[
                            y - 0.18,
                            y + 0.18,
                            y + 0.18,
                            y - 0.18,
                        ],
                        line=dict(width=1, color="#374151"),
                        fillcolor="#374151",
                        opacity=0.8,
                        fill="toself",
                        mode="lines",
                        name=row[roles["feature"]],
                        text=row[roles["ensembl"]],
                        showlegend=False,
                    ),
                )

            elif ftype == "exon":
                y = track_y.get("Transcript")
                if y is None:
                    continue

                head_x = 0  # (stop - start) / 4
                head_y = 0.2
                off_y = 0.2
                new_fig.add_trace(
                    go.Scatter(
                        x=[
                            start,
                            start + head_x,
                            start + head_x,
                            stop,
                            stop,
                            start + head_x,
                            start + head_x,
                            start,
                        ],
                        y=[
                            y,
                            y + head_y,
                            y + off_y,
                            y + off_y,
                            y - off_y,
                            y - off_y,
                            y - head_y,
                            y,
                        ],
                        fill="toself",
                        mode="lines",
                        line=dict(color="black", width=0),
                        fillcolor="#1AADC4",
                        name=row[roles["feature"]],
                        text=row[roles["ensembl"]],
                        showlegend=False,
                    ),
                )
            elif ftype == "intron":
                y = track_y.get("Transcript")
                if y is None:
                    continue

                new_fig.add_trace(
                    go.Scatter(
                        x=[start, start, stop, stop],
                        y=[y + 0.01, y - 0.01, y - 0.01, y + 0.01],
                        mode="lines",
                        fill="toself",
                        fillcolor="#6b7280",
                        line=dict(width=2, color="#6b7280"),
                        hoveron="points+fills",
                        hoverinfo="name",
                        name=row[roles["feature"]],
                        showlegend=False,
                    )
                )

        # batch sgRNA rendering: one filled-polygon trace + one invisible hover trace ---
        sgrna_df = work[work[roles["track"]] == "sgRNA"].copy()
        if not sgrna_df.empty and "sgRNA" in track_y:
            offset_col = roles.get("track_offset")
            xs_combined: list = []
            ys_combined: list = []
            hover_xs: list = []
            hover_ys: list = []
            valid_idx: list = []

            for idx, row in sgrna_df.iterrows():
                try:
                    s = int(row[roles["start"]])
                    e = int(row[roles["stop"]])
                except Exception:
                    continue
                offset = (
                    float(row[offset_col])
                    if offset_col and offset_col in row.index
                    else 0.0
                )
                y_sgrna = track_y["sgRNA"] - offset * state.track_spacing
                head_x, head_y = state.track_head_x * (e - s), state.track_head_y * 1.0
                # Arrow-head polygon; None separates polygons within one trace
                xs_combined += [s, s + head_x, e, e, s + head_x, s, None]
                ys_combined += [
                    y_sgrna,
                    y_sgrna + head_y,
                    y_sgrna + head_y,
                    y_sgrna - head_y,
                    y_sgrna - head_y,
                    y_sgrna,
                    None,
                ]
                hover_xs.append((s + e) / 2.0)
                hover_ys.append(y_sgrna)
                valid_idx.append(idx)

            # determine per-sgRNA fill colors from the metric column (if configured)
            metric_col = getattr(state, "track_metric", "")
            colorscale = getattr(state, "track_metric_colorscale", "RdYlGn")
            use_metric = bool(
                metric_col and metric_col in sgrna_df.columns and valid_idx
            )

            if use_metric:
                raw = sgrna_df.loc[valid_idx, metric_col].to_numpy(dtype=float)
                vmin, vmax = raw.min(), raw.max()
                span = vmax - vmin if vmax != vmin else 1.0
                normalized = ((raw - vmin) / span).tolist()
                fill_colors = plotly.colors.sample_colorscale(colorscale, normalized)
            else:
                fill_colors = ["#facc15"] * len(valid_idx)

            # One filled-polygon trace per sgRNA (needed for individual fill colors)
            # When no metric is set all are the same color so it's still visually one group.
            for poly_idx, (s_idx, idx) in enumerate(
                zip(range(0, len(xs_combined), 7), valid_idx)
            ):
                new_fig.add_trace(
                    go.Scatter(
                        x=xs_combined[s_idx : s_idx + 7],
                        y=ys_combined[s_idx : s_idx + 7],
                        fill="toself",
                        mode="lines",
                        line=dict(color="black", width=1),
                        fillcolor=fill_colors[poly_idx],
                        hoverinfo="skip",
                        name="sgRNA",
                        showlegend=False,  # (poly_idx == 0),
                        legendgroup="sgrna",
                        opacity=state.track_opacity,
                    )
                )

            # Invisible hover markers — one point per sgRNA at its midpoint
            hovertemplate, customdata = make_hover_template(layer, sgrna_df, valid_idx)
            new_fig.add_trace(
                go.Scatter(
                    x=hover_xs,
                    y=hover_ys,
                    mode="markers",
                    marker=dict(size=16, opacity=0, color="rgba(0,0,0,0)"),
                    customdata=customdata,
                    hovertemplate=hovertemplate,
                    showlegend=False,
                    name="sgRNA",
                    legendgroup="sgrna",
                )
            )

            # Dummy colorbar trace so Plotly shows the scale in the figure
            if use_metric:
                raw_all = sgrna_df.loc[valid_idx, metric_col].to_numpy(dtype=float)
                new_fig.add_trace(
                    go.Scatter(
                        x=[None],
                        y=[None],
                        mode="markers",
                        marker=dict(
                            colorscale=colorscale,
                            cmin=float(raw_all.min()),
                            cmax=float(raw_all.max()),
                            color=[float(raw_all.min())],
                            colorbar=dict(title=metric_col, thickness=12),
                            showscale=True,
                            size=0,
                        ),
                        legendgroup="sgrna",
                        showlegend=False,
                        hoverinfo="skip",
                        name=f"{metric_col} scale",
                    )
                )

        new_fig.update_yaxes(
            tickmode="array",
            tickvals=[track_y[name] for name in track_order],
            ticktext=track_order,
            range=[0.5, n_tracks + 0.5],
            showgrid=False,
            zeroline=False,
        )
        xrange = resolve_coordinates(state.track_coordinates, work) or [xmin, xmax]
        new_fig.update_xaxes(
            range=xrange,
            showgrid=True,
            gridcolor="#e5e7eb",
            title="Genomic Position",
        )
        new_fig.update_layout(
            template="plotly_white",
            height=max(380, 120 + 110 * n_tracks),
            margin=dict(l=70, r=30, t=40, b=60),
            coloraxis_colorbar=dict(x=1.08, y=0.5, len=0.75),  # weiter nach rechts
        )

        new_fig = processor.apply(new_fig, state, self.trace_capabilities)
        return _merge(fig, new_fig)


# ─────────────────────────────────────────────────────────────────────────────
# HorizontalBarPlot
# ─────────────────────────────────────────────────────────────────────────────


class HBarPlotState(PlotState):
    sidebar_group = "Layer"

    bar_color = param.Color(default="#1f77b4", label="Bar Color")
    opacity = param.Number(default=1.0, bounds=(0.0, 1.0), label="Opacity")
    top = param.Integer(default=10, bounds=(-100, 100), label="Top N")
    # barmode = param.Selector(
    #     default="group",
    #     objects=["relative", "group", "overlay"],
    #     label="Mode",
    # )
    # aggregation = param.Selector(
    #     default="median",
    #     objects=["sum", "mean", "median"],
    #     label="Aggregation",
    # )


class HorizontalBarPlot(BasePlot):
    ROLE_SPECS = {
        "x": RoleSpec(name="x", dtype="continuous", description="Value axis"),
        "y": RoleSpec(name="y", dtype="categorical", description="Category axis"),
        "error_x": RoleSpec(
            name="error_y",
            dtype="continuous",
            required=False,
            description="Error bar values",
        ),
        "color": RoleSpec(
            name="color",
            dtype="categorical",
            required=False,
            description="Color grouping variable",
        ),
        "group": RoleSpec(
            name="group",
            dtype="categorical",
            required=False,
            description="Grouping variable for stacked bars",
        ),
    }
    state_cls = HBarPlotState
    trace_capabilities: set[str] = set()  # bars don't use marker.* the same way

    def render(
        self,
        df: pd.DataFrame,
        layer: PlotLayer,
        state: HBarPlotState,
        processor: TraceProcessor,
        fig=None,
    ) -> go.Figure:
        roles = layer.roles
        kwargs = {}
        if "color" in roles:
            kwargs["color"] = roles["color"]
            color_map = layer.color_maps.get("color")
            if color_map:
                kwargs["color_discrete_map"] = color_map
        if layer.hover:
            kwargs["hover_data"] = layer.hover
        new_fig = px.bar(
            df,
            x=roles["x"],
            y=roles["y"],
            orientation="h",
            error_x=roles.get("error_x"),
            **kwargs,
        )
        if "color" not in roles:
            new_fig.update_traces(marker_color=state.bar_color)
        new_fig.update_traces(opacity=state.opacity)
        new_fig = processor.apply(new_fig, state, self.trace_capabilities)
        return _merge(fig, new_fig)


################################################################################
# Volcano
################################################################################


class VolcanoPlotState(PlotState):
    """Layout knobs for VolcanoPlot. NOT data roles — those are in PlotLayer."""

    sidebar_group = "Layer"
    x_score = param.Selector(
        default="log2FC",
        objects=["log2FC"],
        label="X-Axis score",
    )

    y_score = param.Selector(
        default="-log10 P-value",
        objects=["FDR", "P-value", "-log10 FDR", "-log10 P-value"],
        label="Y-Axis score",
    )
    fdr_threshold = param.Number(default=0.05, bounds=(0.0, 1.0), label="FDR Threshold")
    logfc_threshold = param.Number(
        default=1.0, bounds=(0.0, 10.0), label="Log2 Fold Change Threshold"
    )
    highlight_significant = param.Boolean(
        default=True, label="Highlight Significant Points"
    )
    marker_size = param.Integer(default=4, bounds=(1, 50), label="Marker Size")
    marker_linewidth = param.Number(
        default=0, bounds=(0.0, 10.0), label="Marker Line Width"
    )
    marker_edgecolor = param.Color(default="#000000", label="Marker Line Color")
    marker_opacity = param.Number(default=1.0, bounds=(0.0, 1.0), label="Opacity")
    threshold_linewidth = param.Number(
        default=1.0, bounds=(0.0, 10.0), label="Threshold Line Width"
    )
    threshold_line = param.Selector(
        default="dash",
        objects=["solid", "dot", "dash", "longdash", "dashdot", "longdashdot"],
        label="Threshold Line Style",
    )
    show_threshold_lines = param.Boolean(default=True, label="Show Threshold Lines")
    positive_color = param.Color(default="#d62728", label="Upregulated Color")
    negative_color = param.Color(default="#1f77b4", label="Downregulated Color")
    neutral_color = param.Color(default="#7f7f7f", label="Neutral Color")
    section = None


class VolcanoPlot(BasePlot):

    ROLE_SPECS = {
        "x": RoleSpec(name="x", dtype="continuous", description="Log2 fold change"),
        "y": RoleSpec(name="y", dtype="continuous", description="-log10 p-value"),
        "color": RoleSpec(
            name="color",
            dtype="categorical",
            required=False,
            description="Optional grouping variable",
        ),
    }

    state_cls = VolcanoPlotState

    trace_capabilities = {
        "marker.size",
        "marker.line.width",
        "marker.line.color",
        "marker.opacity",
        "marker.color",
    }

    def render(
        self,
        df: pd.DataFrame,
        layer: PlotLayer,
        state: VolcanoPlotState,
        processor: TraceProcessor,
        fig=None,
    ) -> go.Figure:
        roles = layer.roles
        kwargs = {}

        x = self.get_param_from_roles(roles["x"], state)
        y = self.get_param_from_roles(roles["y"], state)
        # color scheme
        if "color" in roles:
            color = self.get_param_from_roles(roles["color"], state)
            kwargs["color"] = color
            # color_map = layer.color_maps.get("color")
            color_map = {
                "up": state.positive_color,
                "down": state.negative_color,
                "not significant": state.neutral_color,
            }
            kwargs["color_discrete_map"] = color_map
        # add hover
        if layer.hover:
            kwargs["hover_data"] = layer.hover
        # scatter
        new_fig = px.scatter(df, x=x, y=y, title="Volcano", **kwargs)

        # threshold lines
        if state.show_threshold_lines:
            new_fig.add_shape(
                type="line",
                x0=-state.logfc_threshold,
                x1=-state.logfc_threshold,
                y0=0,
                y1=df[y].max(),
                line={"dash": state.threshold_line, "width": state.threshold_linewidth},
            )
            new_fig.add_shape(
                type="line",
                x0=state.logfc_threshold,
                x1=state.logfc_threshold,
                y1=df[y].max(),
                y0=0,
                line={"dash": state.threshold_line, "width": state.threshold_linewidth},
            )

            new_fig.add_shape(
                type="line",
                x0=df[x].min(),
                x1=df[x].max(),
                y0=-np.log10(state.fdr_threshold),
                y1=-np.log10(state.fdr_threshold),
                line={"dash": state.threshold_line, "width": state.threshold_linewidth},
            )
        # marker styling
        new_fig.update_traces(
            marker={
                "size": state.marker_size,
                "line": {
                    "width": state.marker_linewidth,
                    "color": state.marker_edgecolor,
                },
                "opacity": state.marker_opacity,
            }
        )

        # processor hooks
        new_fig = processor.apply(new_fig, state, self.trace_capabilities)

        if state.highlight_significant:  # this must be after the general trace apply
            new_fig.update_traces(
                selector=lambda t: t.name in {"up", "down"},
                marker={
                    "size": state.marker_size * 1.5,
                    "opacity": 1.0,
                },
            )

        return _merge(fig, new_fig)


################################################################################
# Heatmap
################################################################################


class HeatmapPlotState(PlotState):
    """Layout knobs for HeatmapPlot. NOT data roles — those are in PlotLayer."""

    sidebar_group = "Layer"

    gene_filter = param.String(
        default=None,
        label="Gene Filter",
    )
    sample_filter = param.String(
        default=None,
        label="Sample Filter",
    )
    logfc_threshold = param.Number(
        default=1.0,
        bounds=(0.0, 10.0),
        label="Log2 Fold Change Threshold",
    )
    FDR_threshold = param.Number(
        default=0.05,
        bounds=(0.0, 1.0),
        label="FDR Threshold",
    )
    direction_filter = param.Selector(
        default=None,
        objects=["not significant", "up", "down", None],
        label="Comparison Filter",
    )
    colorscale = param.Selector(
        default="RdBu",
        objects=[
            "Viridis",
            "Plasma",
            "Inferno",
            "Magma",
            "Cividis",
            "RdBu",
            "Blues",
            "Reds",
        ],
        label="Color Scale",
    )

    zmin = param.Number(default=None, allow_None=True, label="Z Min")
    zmax = param.Number(default=None, allow_None=True, label="Z Max")
    show_values = param.Boolean(default=False, label="Show Values in Cells")

    x_tickangle = param.Integer(default=0, bounds=(0, 90), label="X Tick Angle")
    section = None


class HeatmapPlot(BasePlot):

    ROLE_SPECS = {
        "x": RoleSpec(
            name="x", dtype="categorical", description="Columns (samples/features)"
        ),
        "y": RoleSpec(
            name="y", dtype="categorical", description="Rows (genes/features)"
        ),
        "z": RoleSpec(name="z", dtype="continuous", description="Matrix values"),
    }

    state_cls = HeatmapPlotState

    trace_capabilities = {
        "colorscale",
        "zmin",
        "zmax",
        "showscale",
    }

    def render(
        self,
        df: pd.DataFrame,
        layer: PlotLayer,
        state: HeatmapPlotState,
        processor: TraceProcessor,
        fig=None,
    ) -> go.Figure:

        roles = layer.roles

        # --- reshape to matrix form ---
        index = self.get_param_from_roles(roles["y"], state)
        columns = self.get_param_from_roles(roles["x"], state)
        values = self.get_param_from_roles(roles["z"], state)
        pivot = df.pivot(
            index=index,
            columns=columns,
            values=values,
        )

        z = pivot.values
        x = pivot.columns.tolist()
        y = pivot.index.tolist()

        # --- heatmap trace ---
        new_fig = go.Figure(
            data=go.Heatmap(
                z=z,
                x=x,
                y=y,
                colorscale=state.colorscale,
                zmin=state.zmin,
                zmax=state.zmax,
                colorbar={"title": "values"},
                title="Heatmap",
            )
        )

        # --- optional annotations ---
        if state.show_values and z.size <= 2500:
            annotations = []
            for i, yi in enumerate(y):
                for j, xj in enumerate(x):
                    val = z[i][j]
                    if pd.notnull(val):
                        annotations.append(
                            {
                                "x": xj,
                                "y": yi,
                                "text": str(round(val, 2)),
                                "showarrow": False,
                                "font": {"color": "black", "size": 10},
                            }
                        )

            new_fig.update_layout(annotations=annotations)

        # --- layout tweaks ---
        new_fig.update_layout(
            xaxis={"tickangle": state.x_tickangle},
        )

        # --- processor hooks ---
        new_fig = processor.apply(new_fig, state, self.trace_capabilities)

        return _merge(fig, new_fig)


################################################################################
# Dotplot State
################################################################################


class GSEADotplotState(PlotState):

    sidebar_group = "Layer"

    top_x = param.Integer(default=20, bounds=(1, None), label="Top X Terms")
    fdr_threshold = param.Number(
        default=1,
        bounds=(0.0, 1.0),
        label="FDR Threshold",
    )
    geneset_filter = param.String(
        default="",
        label="Filter Gene Sets",
        doc="Filter gene sets by name. For multiple sets provide a comma-separated list.",
    )
    nes_range = param.Range(
        default=(-5.0, 5.0),
        bounds=(-5.0, 5.0),
        doc="Range of normalized enrichment scores (NES) to display",
        label="NES threshold",
    )
    color = param.Selector(
        default="NES",
        objects=["NES", "ES", "FDR q-val", "FWER p-val"],
        label="Color by",
    )
    size = param.Selector(
        default="Leading Edge %",
        objects=["Leading Edge %", "Size", "NES", "ES", "FDR q-val", "FWER p-val"],
        label="Size by",
    )
    sort_by = param.Selector(
        default="NES",
        objects=["FDR q-val", "NES", "ES", "NOM p-val", "FWER p-val"],
        doc="Sort by which column in the GSEA results table.",
        label="Sort By",
    )
    sort_by_ascending = param.Boolean(default=True, label="Sort Ascending")
    color_scale = param.Selector(
        default="Viridis",
        objects=["Viridis", "Cividis", "Plasma", "Inferno", "Magma", "RdBu"],
        label="Color Scale",
    )
    show_colorbar = param.Boolean(default=True, label="Show Colorbar")
    size_max = param.Number(
        default=15,
        bounds=(1, 50),
        label="Size Max",
    )
    y_tickangle = param.Integer(default=0, bounds=(-75, 75), label="Y Tick Angle")

    section = None


################################################################################
# Dotplot Plot
################################################################################


class GSEADotplotPlot(BasePlot):

    ROLE_SPECS = {
        "x": RoleSpec(name="x", dtype="categorical", description="Gene set or term"),
        "y": RoleSpec(
            name="y", dtype="categorical", description="Category (e.g. condition)"
        ),
        "size": RoleSpec(
            name="size", dtype="continuous", description="Gene ratio or count"
        ),
        "color": RoleSpec(
            name="color", dtype="continuous", description="Adjusted p-value or NES"
        ),
    }

    state_cls = GSEADotplotState

    trace_capabilities = {
        "marker.size",
        "marker.color",
        "marker.colorscale",
        "marker.showscale",
    }

    def render(self, df, layer, state, processor, fig=None):

        roles = layer.roles
        x = self.get_param_from_roles(roles["x"], state)
        y = self.get_param_from_roles(roles["y"], state)
        size = self.get_param_from_roles(roles["size"], state)
        color = self.get_param_from_roles(roles["color"], state)
        kwargs = {}
        if layer.hover:
            kwargs["hover_data"] = layer.hover
        new_fig = px.scatter(
            df,
            x=x,
            y=y,
            size=size,
            color=color,
            color_continuous_scale=state.color_scale.lower(),
            range_color=(df[color].min(), df[color].max()),
            size_max=state.size_max,
            title="GSEA",
            **kwargs,
        )
        if not state.show_colorbar:
            new_fig.update_coloraxes(showscale=False)
        new_fig.update_traces(
            marker={
                # size=df["_marker_size"],
                "showscale": state.show_colorbar,
            },
            cliponaxis=False,
        )
        new_fig.update_layout(
            uniformtext_minsize=8,
            uniformtext_mode="hide",
            showlegend=state.show_colorbar,
            yaxis={"tickangle": state.y_tickangle},
        )

        new_fig = processor.apply(new_fig, state, self.trace_capabilities)

        return _merge(fig, new_fig)


################################################################################
# Barplot State
################################################################################


class GSEABarplotState(PlotState):

    sidebar_group = "Layer"

    top_x = param.Integer(default=20, bounds=(1, None), label="Top X Terms")
    fdr_threshold = param.Number(
        default=1,
        bounds=(0.0, 1.0),
        label="FDR Threshold",
    )
    geneset_filter = param.String(
        default="",
        label="Filter Gene Sets",
        doc="Filter gene sets by name. For multiple sets provide a comma-separated list.",
    )
    nes_range = param.Range(
        default=(-5.0, 5.0),
        bounds=(-5.0, 5.0),
        doc="Range of normalized enrichment scores (NES) to display",
        label="NES threshold",
    )
    score = param.Selector(
        default="NES",
        objects=["NES", "ES", "FDR q-val", "FWER p-val"],
        label="Score to display",
    )

    bar_color = param.Selector(
        default="Leading Edge %",
        objects=["Leading Edge %", "Size", "NES", "ES", "FDR q-val", "FWER p-val"],
        label="Color by",
    )
    orientation = param.Selector(
        default="h",
        objects=["h", "v"],
        label="Orientation",
    )
    sort_by = param.Selector(
        default="NES",
        objects=["FDR q-val", "NES", "ES", "NOM p-val", "FWER p-val"],
        doc="Sort by which column in the GSEA results table.",
        label="Sort By",
    )
    sort_by_ascending = param.Boolean(default=True, label="Sort Ascending")
    y_tickangle = param.Integer(default=0, bounds=(-75, 75), label="Y Tick Angle")
    show_values = param.Boolean(default=True, label="Show Values")
    color_scale = param.Selector(
        default="Viridis",
        objects=[
            "Viridis",
            "Balance",
            "Curl",
            "Cividis",
            "Plasma",
            "Inferno",
            "Magma",
            "RdBu",
        ],
        label="Color Scale",
    )
    show_colorbar = param.Boolean(default=True, label="Show Colorbar")
    section = None


################################################################################
# Barplot Plot
################################################################################


class GSEABarplotPlot(BasePlot):

    ROLE_SPECS = {
        "x": RoleSpec(name="x", dtype="categorical", description="Gene set"),
        "y": RoleSpec(
            name="y", dtype="continuous", description="NES or enrichment score"
        ),
        "color": RoleSpec(name="color", dtype="categorical", required=False),
    }

    state_cls = GSEABarplotState

    trace_capabilities = {
        "marker.color",
        "marker.line.width",
        "marker.line.color",
    }

    def render(self, df, layer, state, processor, fig=None):

        roles = layer.roles

        orientation = state.orientation
        x = self.get_param_from_roles(roles["x"], state)
        y = self.get_param_from_roles(roles["y"], state)
        color = self.get_param_from_roles(roles["color"], state)
        kwargs = {}
        if layer.hover:
            kwargs["hover_data"] = layer.hover
        if state.show_values:
            kwargs["text"] = x

        color_continuous_midpoint = None
        if color in ["FDR q-val", "FWER p-val"]:
            color_continuous_midpoint = state.fdr_threshold

        new_fig = px.bar(
            df,
            x=x if orientation == "h" else y,
            y=y if orientation == "h" else x,
            orientation=orientation,
            color=color,
            color_continuous_scale=state.color_scale.lower(),
            color_continuous_midpoint=color_continuous_midpoint,
            title="GSEA Bar",
            **kwargs,
        )
        new_fig = processor.apply(new_fig, state, self.trace_capabilities)
        if not state.show_colorbar:
            new_fig.update_coloraxes(showscale=False)
        if state.show_values:
            new_fig.update_traces(
                texttemplate="%{text:.2f}",
                textposition="inside",
                cliponaxis=False,
            )
        show_text = "show" if state.show_values else "hide"
        new_fig.update_layout(
            showlegend=False,  # state.show_colorbar,
            yaxis={"tickangle": state.y_tickangle},
        )
        new_fig.update_yaxes(
            tickmode="array",
            tickvals=df[y],
            ticktext=[s[: min(len(s), 40)] for s in df[y]],
        )

        return _merge(fig, new_fig)


################################################################################
# GSEA Plot State
################################################################################


class GSEAGseaPlotState(PlotState):

    sidebar_group = "Layer"

    enrichment_color = param.Color(default="#d62728", label="Enrichment Color")

    hit_color = param.Color(default="#000000", label="Hit Color")

    line_width = param.Number(default=2.0, bounds=(0.5, 10.0), label="Line Width")

    show_hits = param.Boolean(default=True, label="Show Gene Hits")

    section = None


################################################################################
# GSEA Plot
################################################################################


class GSEAGseaPlot(BasePlot):

    ROLE_SPECS = {
        "rank": RoleSpec(
            name="rank", dtype="continuous", description="Ranked gene list"
        ),
        "es": RoleSpec(
            name="es", dtype="continuous", description="Enrichment score curve"
        ),
        "hits": RoleSpec(
            name="hits", dtype="continuous", required=False, description="Hit positions"
        ),
    }

    state_cls = GSEAGseaPlotState

    trace_capabilities = {
        "line.color",
        "line.width",
    }

    def render(self, df, layer, state, processor, fig=None):

        roles = layer.roles

        new_fig = go.Figure()

        # enrichment curve
        new_fig.add_trace(
            go.Scatter(
                x=df.index,
                y=df[roles["es"]],
                mode="lines",
                line=dict(
                    color=state.enrichment_color,
                    width=state.line_width,
                ),
                name="Enrichment Score",
            )
        )

        # hits (rug plot style)
        if state.show_hits and "hits" in roles:
            hit_positions = df[roles["hits"]].dropna().values

            new_fig.add_trace(
                go.Scatter(
                    x=hit_positions,
                    y=np.zeros_like(hit_positions),
                    mode="markers",
                    marker=dict(
                        color=state.hit_color,
                        size=4,
                        symbol="line-ns",
                    ),
                    name="Hits",
                )
            )

        new_fig = processor.apply(new_fig, state, self.trace_capabilities)

        return _merge(fig, new_fig)


################################################################################
# PCA Plot State
################################################################################


class PCAPlotState(PlotState):

    sidebar_group = "Layer"

    marker_size = param.Integer(default=8, bounds=(1, 50), label="Marker Size")
    marker_opacity = param.Number(default=0.9, bounds=(0.0, 1.0), label="Opacity")
    show_axes_zero = param.Boolean(default=True, label="Show Zero Axes")

    # color_by_explained_variance = param.Boolean(
    #     default=False, label="Color by Explained Variance"
    # )
    pca_1 = param.Selector(
        default="PC1",
        objects=["PC1", "PC2", "PC3", "PC4", "PC5"],
        doc="Principal component to use for the X-axis.",
        label="X-axis Component",
    )
    pca_2 = param.Selector(
        default="PC2",
        objects=["PC1", "PC2", "PC3", "PC4", "PC5"],
        doc="Principal component to use for the Y-axis.",
        label="Y-axis Component",
    )
    section = None


################################################################################
# PCA Plot
################################################################################


class PCAPlot(BasePlot):

    ROLE_SPECS = {
        "x": RoleSpec(name="x", dtype="continuous", description="PC1 scores"),
        "y": RoleSpec(name="y", dtype="continuous", description="PC2 scores"),
        "color": RoleSpec(
            name="color",
            dtype="categorical",
            required=False,
            description="Optional group label (e.g. condition, batch)",
        ),
        "size": RoleSpec(
            name="size",
            dtype="continuous",
            required=False,
            description="Optional magnitude (e.g. sample weight)",
        ),
    }

    state_cls = PCAPlotState

    trace_capabilities = {
        "marker.size",
        "marker.color",
        "marker.opacity",
        "marker.line.width",
        "marker.line.color",
    }

    def render(
        self,
        df: pd.DataFrame,
        layer,
        state: PCAPlotState,
        processor,
        fig=None,
    ) -> go.Figure:
        roles = layer.roles
        kwargs = {}
        x = self.get_param_from_roles(roles["x"], state)
        y = self.get_param_from_roles(roles["y"], state)
        # grouping
        if "color" in roles:
            color = self.get_param_from_roles(roles["color"], state)
            kwargs["color"] = color

        # optional size
        if "size" in roles:
            size = self.get_param_from_roles(roles["size"], state)
            kwargs["size"] = size

        if layer.hover:
            kwargs["hover_data"] = layer.hover

        # base scatter
        new_fig = px.scatter(df, x=x, y=y, title="PCA", **kwargs)

        # styling
        new_fig.update_traces(
            marker=dict(
                size=state.marker_size,
                opacity=state.marker_opacity,
            )
        )

        # zero axes (classic PCA crosshair)
        if state.show_axes_zero:
            x_min, x_max = df[x].min(), df[x].max()
            y_min, y_max = df[y].min(), df[y].max()

            new_fig.add_shape(
                type="line",
                x0=0,
                x1=0,
                y0=y_min,
                y1=y_max,
                line=dict(width=1, dash="dot"),
            )

            new_fig.add_shape(
                type="line",
                x0=x_min,
                x1=x_max,
                y0=0,
                y1=0,
                line=dict(width=1, dash="dot"),
            )

        # processor hook (your global styling system)
        new_fig = processor.apply(new_fig, state, self.trace_capabilities)

        return _merge(fig, new_fig)


################################################################################
# Scree
################################################################################
class ScreePlotState(PlotState):
    show_cumulative = param.Boolean(default=True)
    n_components = param.Integer(default=5, bounds=(1, None))
    bar_color = param.Color(default="#4C78A8")
    cumulative_color = param.Color(default="#E45756")


class ScreePlot(BasePlot):
    ROLE_SPECS = {
        "component": RoleSpec(
            name="component",
            dtype="categorical",
            description="PC label (e.g. PC1, PC2, ...)",
        ),
        "variance": RoleSpec(
            name="variance",
            dtype="continuous",
            description="Explained variance ratio per PC",
        ),
    }
    state_cls = ScreePlotState
    trace_capabilities = {
        "marker.color",
        "marker.line.width",
        "marker.line.color",
        "marker.opacity",
    }

    def render(
        self,
        df: pd.DataFrame,
        layer,
        state: ScreePlotState,
        processor,
        fig=None,
    ) -> go.Figure:
        roles = layer.roles
        dff = df.head(state.n_components)

        new_fig = go.Figure()
        new_fig.add_bar(
            x=dff[roles["component"]],
            y=dff[roles["variance"]],
            marker=dict(color=state.bar_color),
            name="Explained variance",
        )

        if state.show_cumulative:
            cum = dff[roles["variance"]].cumsum()
            new_fig.add_scatter(
                x=dff[roles["component"]],
                y=cum,
                mode="lines+markers",
                line=dict(color=state.cumulative_color),
                name="Cumulative",
                yaxis="y2",
            )
            new_fig.update_layout(
                yaxis2=dict(
                    overlaying="y",
                    side="right",
                    range=[0, 1.05],
                    title="Cumulative variance",
                )
            )

        new_fig = processor.apply(new_fig, state, self.trace_capabilities)
        return _merge(fig, new_fig)


################################################################################
# convenience
################################################################################


def _merge(fig: go.Figure | None, new_fig: go.Figure) -> go.Figure:
    if fig is None:
        return new_fig
    for trace in new_fig.data:
        fig.add_trace(trace)
    return fig


def build_pca_overview(
    pca_scores_df: pd.DataFrame, variance_df: pd.DataFrame
) -> pn.viewable.Viewable:
    scatter_plot = PanelPlot(
        plot_cls=PCAPlot,
        df=pca_scores_df,
        layer=Layer(roles={"x": "PC1", "y": "PC2", "color": "condition"}),
    )
    scree_plot = PanelPlot(
        plot_cls=ScreePlot,
        df=variance_df,
        layer=Layer(roles={"component": "component", "variance": "variance"}),
    )
    return pn.Row(scatter_plot.panel(), scree_plot.panel())


# ─────────────────────────────────────────────────────────────────────────────
# Cross-cutting shared state — independent of which plots are combined
# ─────────────────────────────────────────────────────────────────────────────


class TitleState(PlotState):
    sidebar_group = "Appearance"

    appearance_title = param.String(default="", label="Title")
    appearance_xlabel = param.String(default=None, allow_None=True, label="X Label")
    appearance_ylabel = param.String(default=None, allow_None=True, label="Y Label")


class ExportState(PlotState):
    sidebar_group = "Export"

    export_folder = param.String(default="results", label="Export Folder")
    export_filename = param.String(default="figure", label="File Name")

    # actions
    export_svg = param.Action(action_export_figure("svg"), label="Export SVG")
    export_png = param.Action(action_export_figure("png"), label="Export PNG")
    export_pdf = param.Action(action_export_figure("pdf"), label="Export PDF")

    export_tables = param.Action(action_export_data(format="tsv"), label="Export TSV")
    # export_bed = param.Action(lambda self: None, label="Export BED")
    # export_bigwig = param.Action(lambda self: None, label="Export BigWig")
    export_app = param.Action(action_export_app(), label="Export Standalone App")

    # reset_button = param.Action(
    #     reset(
    #         [
    #             f"{sidebar_group}_select_sample",
    #             f"{sidebar_group}_select_mouse",
    #             f"{sidebar_group}_select_tissue",
    #             f"{sidebar_group}_select_trp53_status",
    #         ]
    #     ),
    #     label="Reset Selection",
    # )


# ─────────────────────────────────────────────────────────────────────────────
# Cross-cutting shared state — Data Selection
# ─────────────────────────────────────────────────────────────────────────────


class CellSelectionState(PlotState):
    analysis_cellline = param.Selector(
        default=None,
        objects=[None, "HCT_D4", "HuH7"],
        doc="Cell line to filter the data by. If None, no filtering is applied.",
        label="Cell Line",
    )


class FilterSelectionState(PlotState):
    sidebar_group = "Filters"
    analysis_filter = param.String(
        default="",
        label="Filter",
        doc=r'Filter expression to apply to the data. If empty, no filtering is applied. Filter expressions can be something like: `sgRNA == "GGACCGGCGCACAGAGGAAG"`, `Count > 100` or `Count == 100` or \`sgRNA Number\` == "N0".',
    )


class TreatmentSelectionState(DynamicPlotState):
    sidebar_group = "Filters"
    DYNAMIC_SELECTORS = [
        DynamicSelectorSpec(param_name="analysis_treatment", role="treatment"),
    ]
    # objects=[None] is just a placeholder until populate_dynamic_selectors()
    # overwrites it with the REAL treatments from the data (Rez, DMSO, ...)
    analysis_treatment = param.Selector(
        default=None,
        objects=[None],
        label="Treatment",
        doc="Treatment to filter the data by. If None, no filtering is applied.",
    )


class FDRSelectionState(PlotState):

    analysis_fdr = param.Number(
        default=0.05,
        bounds=(0, 1),
        doc="False discovery rate threshold.",
        label="FDR Threshold",
    )


# ─────────────────────────────────────────────────────────────────────────────
# Registry — plot_type string -> class
# ─────────────────────────────────────────────────────────────────────────────

PLOT_REGISTRY: dict[str, type[BasePlot]] = {
    "rate": RatePlot,
    "hbar": HorizontalBarPlot,
    "genetracks": GenePlot,
    "volcano": VolcanoPlot,
    "heat": HeatmapPlot,
    "gsea": GSEADotplotPlot,
    "gsea_bar": GSEABarplotPlot,
    "gsea_enrich": GSEAGseaPlot,
    "pca": PCAPlot,
    "scree": ScreePlot,
}


def register_plot(name: str, cls: type[BasePlot]) -> None:
    PLOT_REGISTRY[name] = cls


def assign_y(
    df: DataFrame, start: str = "start", stop: str = "stop", y_col: str = "y"
) -> Series:
    """
    Greedy Min-Heap interval scheduling algorithm to assign y-values to intervals in a
    DataFrame.

    Parameters
    ----------
    df : DataFrame
        DataFrame containing the intervals.
    start : str
        Column name for the start of the intervals, by default "start".
    stop : str
        Column name for the end of the intervals, by default "stop".
    y_col : str
        Column name for the y-values to be assigned, by default "y".

    Returns
    -------
    Series
        Series with the y-values for each interval, indicating the
        assigned lane for each interval.
    """
    df = df.sort_values(start)
    heap, free, y, out = [], [], 0, []

    for s, e in zip(df[start], df[stop]):
        while heap and heap[0][0] <= s:
            free.append(heapq.heappop(heap)[1])

        if free:
            lane = free.pop()
        else:
            lane = y
            y += 1

        heapq.heappush(heap, (e, lane))
        out.append(lane)

    return pd.Series(out, index=df.index, name=y_col, dtype=int)


def count_overlaps(
    df: DataFrame,
    start: str = "start",
    stop: str = "stop",
    overlap_col: str = "overlap_count",
) -> Series:
    """
    Count the number of overlapping intervals for each row in a DataFrame.

    Parameters
    ----------
    df : DataFrame
        DataFrame containing the intervals.
    start : str, optional
        Column name for the start of the intervals, by default "start"
    stop : str, optional
        Column name for the end of the intervals, by default "stop"
    overlap_col : str, optional
        Column name for the overlap count, by default "overlap_count"

    Returns
    -------
    Series
        Series with the overlap count for each interval.
    """
    starts = sorted(df[start])
    ends = sorted(df[stop])

    overlap_counts = []

    for _, r in df.iterrows():
        # intervals that started before r.stop
        a = bisect.bisect_right(starts, r[stop])
        # intervals that ended before r.start
        b = bisect.bisect_left(ends, r[start])

        overlap_counts.append(a - b)

    return pd.Series(overlap_counts, index=df.index, name=overlap_col, dtype=int)


def resolve_coordinates(
    coord_str: str,
    df: DataFrame,
    start_col: str = "Start",
    end_col: str = "Stop",
    feature_col: str = "Feature",
    type_col: str = "Type",
) -> tuple[float, float] | None:
    if not coord_str:
        return None
    elif coord_str in df[feature_col].values:
        feature_row = df[df[feature_col] == coord_str].iloc[0]
        start = float(feature_row[start_col])
        end = float(feature_row[end_col])
        return (start, end)
    elif coord_str in df[type_col].values:
        feature_df = df[df[type_col] == coord_str]
        start = (
            feature_df.loc[feature_df[start_col] > 0, [start_col, end_col]].min().min()
        )
        end = (
            feature_df.loc[feature_df[start_col] > 0, [start_col, end_col]].max().max()
        )
        start = float(feature_df[start_col].min())
        end = float(feature_df[end_col].max())
        return (start, end)

    elif "-" in coord_str:
        parts = coord_str.split("-")
        if len(parts) != 2:
            return None
        start_str, end_str = parts
        start = float(start_str.replace(",", "").strip())
        end = float(end_str.replace(",", "").strip())
        return (start, end)
    else:
        raise ValueError(f"Invalid coordinates: {coord_str}")


def make_hover_template(
    layer: PlotLayer, df: DataFrame, valid_idx: list[Any]
) -> tuple[str, np.ndarray | None]:
    hover_cols = [c for c in (layer.hover or []) if c in df.columns]
    customdata = df.loc[valid_idx, hover_cols].to_numpy() if hover_cols else None
    hovertemplate = (
        "<br>".join(f"{c}: %{{customdata[{i}]}}" for i, c in enumerate(hover_cols))
        + "<extra></extra>"
        if hover_cols
        else "Position: %{x:.0f}<extra></extra>"
    )
    return hovertemplate, customdata
