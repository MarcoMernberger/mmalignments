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

import pandas as pd
import param
import plotly.express as px
import plotly.graph_objects as go

from .spec import (
    BasePlot,
    PlotLayer,
    PlotState,
    RoleSpec,
    DynamicPlotState,
    DynamicSelectorSpec,
)

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
        if getattr(combined_state, "appearance_xlabel", None):
            fig.update_xaxes(title=combined_state.appearance_xlabel)
        if getattr(combined_state, "appearance_ylabel", None):
            fig.update_yaxes(title=combined_state.appearance_ylabel)
        if getattr(combined_state, "appearance_title", None):
            fig.update_layout(title=combined_state.appearance_title)
        return fig

    def _update_tick_labels(
        self, fig, combined_state, layers: list[PlotLayer], result
    ) -> go.Figure:
        """Relabel x-ticks using a different column than the one plotted
        (e.g. plot numeric `position`, show `sample` as tick text). Uses
        the FIRST layer's source frame as the reference for tick values."""
        xtick_col = getattr(combined_state, "appearance_xtick", None)
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


def _merge(fig: go.Figure | None, new_fig: go.Figure) -> go.Figure:
    if fig is None:
        return new_fig
    for trace in new_fig.data:
        fig.add_trace(trace)
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Cross-cutting shared state — independent of which plots are combined
# ─────────────────────────────────────────────────────────────────────────────


class TitleState(PlotState):
    sidebar_group = "General"

    appearance_title = param.String(default="", label="Title")
    appearance_xlabel = param.String(default=None, allow_None=True, label="X Label")
    appearance_ylabel = param.String(default=None, allow_None=True, label="Y Label")


class ExportState(PlotState):
    sidebar_group = "Export"

    export_folder = param.String(default="results", label="Export Folder")
    export_filename = param.String(default="figure", label="File Name")
    export_png = param.Action(lambda self: None, label="Export PNG")
    export_svg = param.Action(lambda self: None, label="Export SVG")
    export_tables = param.Action(lambda self: None, label="Export TSV")
    export_app = param.Action(lambda self: None, label="Export Standalone App")


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
}


def register_plot(name: str, cls: type[BasePlot]) -> None:
    PLOT_REGISTRY[name] = cls
