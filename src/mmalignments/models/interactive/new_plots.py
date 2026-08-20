from __future__ import annotations

import plotly.express as px  # type: ignore[import]
import plotly.graph_objects as go  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

from mmalignments.models.interactive.spec import (
    LayerStateView,
)

from .spec import (
    BasePlot,
    PlotLayer,
)


class BasePlot:
    pass


class BarPlot(BasePlot):
    ROLES = ["color", "x", "y", "error", "group", "orientation"]
    state_cls = HBarPlotState
    trace_capabilities: set[str] = set()  # bars don't use marker.* the same way

    def render(
        self,
        df: DataFrame,  # the data
        state: Plotstate,  # defines the params
        layer: PlotLayer,  # the plot layer
        view: LayerStateView,  # resolves roles
        processor: TraceProcessor,
        fig=None,
    ) -> go.Figure:
        roles = view.resolve_roles(self.ROLES)
        kwargs = {}
        x = roles["x"]
        y = roles["y"]
        orientation = view["orientation"]
        color = view[roles["color"]]
        error = view[roles["error"]] if "error" in roles else None
        if "color" in roles:
            kwargs["color"] = color
            color_map = layer.color_maps.get("color")
            if color_map:
                kwargs["color_discrete_map"] = color_map
        if layer.hover:
            kwargs["hover_data"] = layer.hover
        new_fig = px.bar(
            df,
            x=x,
            y=y,
            orientation=orientation,
            error_x=error if orientation == "h" else None,
            error_y=error if orientation == "v" else None,
            **kwargs,
        )
        if "color" not in roles:
            new_fig.update_traces(marker_color=state.bar_color)
        for param in state.param:
            print(param)
        new_fig.update_traces(opacity=view["opacity"])
        new_fig = processor.apply(new_fig, state, self.trace_capabilities)
        return _merge(fig, new_fig)
