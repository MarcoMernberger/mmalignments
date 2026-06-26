import ast
import inspect
import subprocess
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import panel as pn  # type: ignore[import]
import param  # type: ignore[import]
import plotly.express as px  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]
from panel.pane import Plotly  # type: ignore[import]
from panel.viewable import Viewable  # type: ignore[import]
from param import Parameterized  # type: ignore[import]

from mmalignments.services.io import parents  # type: ignore[import]

pn.extension("plotly")


@dataclass
class AnalysisResult:
    raw: DataFrame
    agg: DataFrame


class ColorPickerUI:
    def __init__(self, state):
        self.state = state
        self._widgets = {}

    def build(self, treatments):
        rows = []

        for t in sorted(treatments):
            current = self.state.color_map.get(t, "#cccccc")

            picker = pn.widgets.ColorPicker(
                name=t,
                value=current,
                width=140,
            )

            def _update(event, t=t):
                cm = dict(self.state.color_map)
                cm[t] = event.new
                self.state.color_map = cm

            picker.param.watch(_update, "value")

            self._widgets[t] = picker
            rows.append(pn.Row(pn.pane.Str(t, width=120), picker))

        return pn.Column(*rows)


# ── State ─────────────────────────────────────────────────────────────────────


class BasePlotState(Parameterized):
    """
    All tunable parameters — analysis + layout.

    Contains all tunable parameters that can be interactively accessed in the
    plot and will trigger a re-render. This is a single flat Parameterized class
    that combines all parameters for both analysis and layout.
    For each specific plot, you would create a subclass of PlotState with the
    relevant parameters. The parameters can be grouped logically by prefixing
    their names (e.g., "analysis_" for analysis parameters, "plot_" for layout
    parameters), and the UI can be organised accordingly in the Panel layout.
    The export_app() method will automatically serialise the current state as
    defaults for the app.
    """

    # analysis
    analysis_normalize = param.Boolean(
        default=False, doc="Normalize values to [0,1]"
    )  # noqa: E501
    analysis_window_min = param.Number(default=0.0, doc="Lower bound filter")
    analysis_window_max = param.Number(default=100.0, doc="Upper bound filter")

    # plot layout
    plot_marker_size = param.Integer(default=8, bounds=(1, 40))
    plot_color_by = param.Selector(default="group", objects=["group", "value"])
    plot_title = param.String(default="Scatter")
    plot_legend_fontsize = param.Integer(default=11, bounds=(6, 24))

    # io
    io_output_dir = param.String(default="results")
    io_stem = param.String(default="analysis")


class PlotState(Parameterized):
    """
    All tunable parameters — analysis + layout.

    Contains all tunable parameters that can be interactively accessed in the
    plot and will trigger a re-render. This is a single flat Parameterized class
    that combines all parameters for both analysis and layout.
    For each specific plot, you would create a subclass of PlotState with the
    relevant parameters. The parameters can be grouped logically by prefixing
    their names (e.g., "analysis_" for analysis parameters, "plot_" for layout
    parameters), and the UI can be organised accordingly in the Panel layout.
    The export_app() method will automatically serialise the current state as
    defaults for the app.
    """

    # analysis
    analysis_normalize = param.Boolean(
        default=False, doc="Normalize values to [0,1]"
    )  # noqa: E501
    analysis_window_min = param.Number(default=0.0, doc="Lower bound filter")
    analysis_window_max = param.Number(default=100.0, doc="Upper bound filter")

    # plot layout
    plot_marker_size = param.Integer(default=8, bounds=(1, 40))
    plot_color_by = param.Selector(default="group", objects=["group", "value"])
    plot_title = param.String(default="Scatter")
    plot_legend_fontsize = param.Integer(default=11, bounds=(6, 24))

    # io
    io_output_dir = param.String(default="results")
    io_stem = param.String(default="analysis")


# class EditingRatePlotState(Parameterized):
#     # analysis
#     analysis_metric = param.Selector(
#         default="editing_rate",
#         objects=["editing_rate", "delta_editing_rate"],
#     )

#     analysis_gene = param.String(default="")

#     analysis_show_untreated = param.Boolean(
#         default=True,
#         doc="Include untreated samples",
#     )

#     analysis_plot_type = param.Selector(
#         default="scatter",
#         objects=["scatter", "bar", "violin"],
#     )

#     # appearance
#     plot_title = param.String(default="Editing Rate")

#     plot_title_fontsize = param.Integer(
#         default=18,
#         bounds=(6, 60),
#     )

#     plot_label_fontsize = param.Integer(
#         default=14,
#         bounds=(6, 40),
#     )

#     plot_legend_fontsize = param.Integer(
#         default=12,
#         bounds=(6, 40),
#     )

#     plot_marker_size = param.Integer(
#         default=8,
#         bounds=(1, 50),
#     )

#     plot_line_width = param.Number(
#         default=1.0,
#         bounds=(0.0, 10.0),
#     )

#     plot_facecolor = param.Color(default="#ffffff")
#     plot_edgecolor = param.Color(default="#000000")

#     # treatment colors
#     color_AV_01 = param.Color(default="#1f77b4")
#     color_AV_02 = param.Color(default="#ff7f0e")
#     color_AV_03 = param.Color(default="#2ca02c")
#     color_untreated = param.Color(default="#999999")

#     # export
#     io_output_dir = param.String(default="results")
#     io_stem = param.String(default="editing_rates")

#     io_export_png = param.Action(lambda self: None)
#     io_export_svg = param.Action(lambda self: None)
#     io_export_pdf = param.Action(lambda self: None)


# ── Analysis (pure compute, no UI) ────────────────────────────────────────────


class BaseAnalysis:
    """Override run() to transform raw data given the current state."""

    def run(self, data: DataFrame, state: Parameterized) -> DataFrame:
        raise NotImplementedError


class InteractiveAnalysis(BaseAnalysis):
    """
    Stateless transformer: data + state → derived DataFrame.

    This takes care of the input data, transforms it into a tidy format, and
    applies any filters or normalisation based on the state. The resulting
    DataFrame must have columns 'x', 'y', and any others needed for plotting
    (e.g., 'group' for color).
    """

    def run(self, data: DataFrame, state: Parameterized) -> DataFrame:
        """
        Perform analysis and return a tidy DataFrame ready for plotting. This
        takes care of all data transformations, filtering, and normalisation
        based on the current state. It must contain all processing that should
        be interactively applied and reflected in the plot. Any non-interactive
        processing should be done outside of this method.
        Parameters
        ----------
        data : DataFrame
            Input data to be analyzed.
        state : PlotState
            Current state containing all tunable parameters.

        Returns
        -------
        DataFrame
            A tidy DataFrame with columns 'x', 'y', and any others needed for
            plotting (e.g., 'group' for color).
        """
        data["transformed_value"] = data[
            data["value"].between(
                state.analysis_window_min, state.analysis_window_max
            )  # noqa: E501
        ]
        if state.analysis_normalize:
            v = data["value"]
            data["value"] = (v - v.min()) / (v.max() - v.min())
        return data


# ── Plot (pure render, no UI) ─────────────────────────────────────────────────


class BasePlot:
    """Override render() to return a Plotly figure."""

    def render(self, df: DataFrame, state: PlotState):
        raise NotImplementedError


class InteractivePlot(BasePlot):
    """
    Stateless renderer: derived DataFrame + state → Plotly figure.

    This contains the actual plot function. For each type of plot we need a new
    subclass of this that implements the render() method. The render() method
    takes the tidy DataFrame produced by the analysis and the current state, and
    returns a Plotly figure. The state can be used to control layout options
    like marker size, color scheme, title, etc.
    """

    def render(self, df: DataFrame, state: PlotState):
        fig = px.scatter(
            df,
            x="x",
            y="y",
            color=state.plot_color_by,
            size_max=state.plot_marker_size,
            title=state.plot_title,
        )
        fig.update_layout(
            legend={"font": {"size": state.plot_legend_fontsize}},
            margin={"l": 40, "r": 20, "t": 40, "b": 40},
        )
        return fig


class EditingRatePlot(BasePlot):

    def _treatment_color_map(self, state):

        return {
            "AV_01": state.color_AV_01,
            "AV_02": state.color_AV_02,
            "AV_03": state.color_AV_03,
            "untreated": state.color_untreated,
        }

    def render(self, df, state):

        dff = df.copy()

        # -------------------------
        # filtering
        # -------------------------

        if state.analysis_gene:
            dff = dff[dff["gene"] == state.analysis_gene]

        if not state.analysis_show_untreated:
            dff = dff[~dff["treatment"].str.contains("untreated", case=False, na=False)]

        if (
            state.analysis_metric == "delta_editing_rate"
            and "untreated" in dff["treatment"].astype(str).str.lower().unique()
        ):
            dff = dff[~dff["treatment"].str.contains("untreated", case=False, na=False)]

        colors = self._treatment_color_map(state)

        # -------------------------
        # scatter
        # -------------------------

        if state.analysis_plot_type == "scatter":

            fig = px.scatter(
                dff,
                x="position",
                y=state.analysis_metric,
                color="treatment",
                symbol="mouse",
                color_discrete_map=colors,
                hover_data=[
                    "mouse",
                    "gene",
                    "sample",
                ],
            )

            fig.update_traces(
                marker=dict(
                    size=state.plot_marker_size,
                    line=dict(
                        width=state.plot_line_width,
                        color=state.plot_edgecolor,
                    ),
                )
            )

        # -------------------------
        # bar
        # -------------------------

        elif state.analysis_plot_type == "bar":

            fig = px.bar(
                dff,
                x="position",
                y=state.analysis_metric,
                color="treatment",
                barmode="group",
                color_discrete_map=colors,
            )

        # -------------------------
        # violin
        # -------------------------

        else:

            fig = px.violin(
                dff,
                x="treatment",
                y=state.analysis_metric,
                color="treatment",
                box=True,
                points="all",
                color_discrete_map=colors,
            )

        # -------------------------
        # layout
        # -------------------------

        fig.update_layout(
            title=state.plot_title,
            title_font_size=state.plot_title_fontsize,
            font_size=state.plot_label_fontsize,
            legend_font_size=state.plot_legend_fontsize,
            paper_bgcolor=state.plot_facecolor,
            plot_bgcolor=state.plot_facecolor,
        )

        return fig


# ── PanelPlot (UI + reactivity) ───────────────────────────────────────────────


class PanelPlot(Parameterized):
    """
    Wires state → analysis → plot, exposes Panel layout.

    This is where everything comes together. PanelPlot takes the input data, an
    instance of InteractiveAnalysis, and an instance of InteractivePlot. It also
    holds the PlotState, which is the single source of truth for all tunable
    parameters. The _figure_pane method is decorated with @param.depends on the
    entire state, so it re-runs whenever any parameter changes. Inside
    _figure_pane, we call the analysis to get the derived DataFrame, then pass
    that to the plot to get a Plotly figure, which is returned as a Panel pane.
    The view() method defines the overall layout of the app, including the
    sidebar with parameter widgets and the main area with the plot.
    Finally, export_app() allows us to export the app as a self-contained Python
    script that can be converted to a standalone HTML file using Panel's CLI
    tool and Pyodide.
    """

    state = param.Parameter()

    def __init__(
        self,
        data: DataFrame,
        analysis: BaseAnalysis,
        plot: BasePlot,
        state: PlotState | None = None,
        **params,
    ):
        super().__init__(**params)
        self.data = data
        self.analysis = analysis
        self.plot = plot
        self.state = state or PlotState()
        self.init_state()
        self._last_fig = None

    # ── reactive core ─────────────────────────────────────────────────────────

    def init_state(self) -> None:
        """change state from data"""
        for state_key, value in self.state.update_from_frame.items():
            if value in self.data.columns:
                unique_values = self.data[value].dropna().unique().tolist()
                param_obj = self.state.param[state_key]
                if isinstance(param_obj, param.ListSelector):
                    param_obj.objects = unique_values
                    param_obj.default = None
            if state_key == "export_png":
                self.state.export_png = lambda s: self._save_png(s)
            elif state_key == "export_svg":
                self.state.export_svg = lambda s: self._save_svg(s)
            elif state_key == "export_pdf":
                self.state.export_pdf = lambda s: self._save_pdf(s)
        self._initialized = True

    @param.depends("state.param")  # re-runs on ANY state change
    def _figure_pane(self) -> Plotly:
        fig = self.__build_figure()
        return pn.pane.Plotly(fig, sizing_mode="stretch_width", height=500)

    def __build_figure(self):
        result = self.analysis.run(
            self.data,
            self.state,
        )
        fig = self.plot.render(
            result,
            self.state,
        )
        self._last_fig = fig
        return fig

    # ── public layout ─────────────────────────────────────────────────────────
    def panel(self) -> Viewable:
        # split params into logical groups
        param_groups = {"Analysis": [], "Appearance": [], "Export": []}
        for group_name, group_params in self.state.ui_groups.items():
            param_groups[group_name].extend(group_params)

        widgets = []
        for group_name, group_params in param_groups.items():
            widgets.append(
                (
                    group_name,
                    pn.Param(
                        self.state,
                        parameters=group_params,
                        show_name=True,
                        name=group_name,
                        # widgets={"export_pdf": pn.widgets.Button},
                    ),
                )
            )

        sidebar = pn.Accordion(
            *widgets,
            active=[0],
        )

        return pn.template.FastListTemplate(
            title="Interactive Analysis",
            sidebar=[
                sidebar,
            ],
            main=[self._figure_pane],
        )

    ############################################################################
    # IO
    ############################################################################
    def _save_png(self, event):
        self.save_figure("png")

    def _save_svg(self, event):
        self.save_figure("svg")

    def _save_pdf(self, event):
        self.save_figure("pdf")

    def save_figure(
        self,
        ext: str,
    ):
        print(ext)
        print("Saving figure...")
        if self._last_fig is not None:
            folder = getattr(
                self.state,
                "export_folder",
                getattr(self.state, "folder", "results"),
            )
            filename = getattr(
                self.state,
                "export_filename",
                getattr(self.state, "name", "figure"),
            )
            outfile = Path(folder) / f"{filename}.{ext}"
            parents(outfile)
            self._last_fig.write_image(f"{outfile}")
            print(f"Figure saved to {outfile}")

            if "export_savefigstatus" in self.state.param:
                self.state.export_savefigstatus = f"Figure saved to {outfile}"

    # ── WASM export ───────────────────────────────────────────────────────────

    def export_app(self, data_path: Path | str | None = None) -> Path:
        """
        Write a self-contained app.py + requirements.txt into output_dir/stem/.
        Run afterwards:
            panel convert <outdir>/<stem>/app.py --to pyodide-worker --out dist/
        """
        outdir = Path(self.state.io_output_dir) / self.state.io_stem
        outdir.mkdir(parents=True, exist_ok=True)

        # serialise current state as defaults for the exported app
        state_defaults: dict[str, Any] = {
            name: getattr(self.state, name)
            for name in self.state.param
            if name != "name"
        }

        # serialise data
        if data_path is None:
            dpath = outdir / "data.parquet"
            print(dpath)
            self.data.to_parquet(dpath, index=False)
        else:
            dpath = Path(data_path)

        # app_code = _render_app_template(
        #     data_path=dpath,
        #     state_cls=type(self.state).__name__,
        #     state_defaults=state_defaults,
        #     analysis_cls=type(self.analysis).__name__,
        #     plot_cls=type(self.plot).__name__,
        # )
        state_type = type(self.state)
        analysis_type = type(self.analysis)
        plot_type = type(self.plot)
        app_code = _render_app_py(
            datapath=str(dpath),
            state_src=_clean_source(state_type),
            analysis_src=_clean_source(analysis_type),
            plot_src=_clean_source(plot_type),
            state_cls_name=state_type.__name__,
            analysis_cls_name=analysis_type.__name__,
            plot_cls_name=plot_type.__name__,
            # data_b64=data_b64,
            param_defaults=state_defaults,
        )

        (outdir / "app.py").write_text(app_code)
        (outdir / "requirements.txt").write_text(
            "panel\nplotly\npandas\nparam\npyarrow\n"
        )

        print(f"App written to {outdir}")
        print(
            f"Convert with:\n  panel convert {outdir}/app.py "
            f"--to pyodide-worker --out {outdir}/dist/"
        )
        return outdir

    def export_html(self, data_path: Path | str | None = None) -> Path:
        """
        Convenience method: export app and convert to stand-alone HTML in one
        step. This requires the panel CLI tool to be installed and on the PATH.
        The result is a single HTML file that can be shared directly and runs
        entirely in the browser via WebAssembly.
        """
        outdir = self.export_app(data_path)
        dist_dir = outdir / "dist"
        dist_dir.mkdir(exist_ok=True)
        cmd = f"panel convert {outdir}/app.py --to pyodide-worker --out {dist_dir}"  # noqa: E501
        print(f"Running:\n  {cmd}")

        subprocess.run(cmd, shell=True, check=True)
        html_path = dist_dir / "app.html"
        print(f"Standalone HTML exported to {html_path}")
        return html_path


# -- widget fiddling (optional) ---
# def build_param_widgets(
#     state: param.Parameterized,
#     group_params: list[str],
#     overrides: dict[str, type] | None = None,
# ):
#     widgets = []
#     print(overrides)
#     for name in group_params:
#         param_obj = state.param[name]

#         widget_cls = overrides.get(name, None) if overrides is not None else None
#         if widget_cls is not None:
#             # custom widget
#             widget = widget_cls(
#                 name=name,
#                 value=getattr(state, name),
#                 options=getattr(param_obj, "objects", None),
#             )
#             widget = bind(widget, state, name)

#         else:

#             # fallback to normal Panel Param behavior
#             widget = pn.Param(state, parameters=[name], show_name=True)

#         print(name, widget_cls, widget)
#         widgets.append(widget)

#     return pn.Column(*widgets)


# def bind(widget, state, attr):
#     def _sync(event):
#         setattr(state, attr, event.new)

#     widget.param.watch(_sync, "value")
#     return widget


# ── App template renderer ─────────────────────────────────────────────────────


def _render_app_template(
    data_path: Path,
    state_cls: str,
    state_defaults: dict,
    analysis_cls: str,
    plot_cls: str,
) -> str:
    """Generates the standalone app.py source as a string."""

    defaults_repr = ",\n    ".join(
        f"{k}={v!r}" for k, v in state_defaults.items()
    )  # noqa: E501

    return f'''"""Auto-generated Panel app — edit carefully."""
import panel as pn
import pandas as pd

pn.extension("plotly")

# ── import your analysis classes (adjust path as needed) ──────────────────────
from my_module import {state_cls}, {analysis_cls}, {plot_cls}, PanelPlot

data = pd.read_parquet("{data_path}")

state    = {state_cls}(
    {defaults_repr}
)
analysis = {analysis_cls}()
plot_obj = {plot_cls}()

app = PanelPlot(data, analysis, plot_obj, state=state)
app.view().servable()
'''


def _clean_source(cls: type) -> str:
    src = inspect.getsource(cls)
    tree = ast.parse(src)

    for node in ast.walk(tree):
        # Klassen-Docstring entfernen
        if isinstance(node, ast.ClassDef):
            if (
                node.body
                and isinstance(node.body[0], ast.Expr)
                and isinstance(node.body[0].value, ast.Constant)
                and isinstance(node.body[0].value.value, str)
            ):
                node.body.pop(0)

    cleaned = ast.unparse(tree)
    return textwrap.dedent(cleaned)


def _render_app_py(
    datapath: str,
    state_src: str,
    analysis_src: str,
    plot_src: str,
    state_cls_name: str,
    analysis_cls_name: str,
    plot_cls_name: str,
    # data_b64: str,
    param_defaults: dict[str, Any],
) -> str:
    defaults_repr = "\n    ".join(f"{k}={v!r}," for k, v in param_defaults.items())

    return f'''\
"""
Auto-generated Panel/Pyodide app.
DO NOT EDIT — regenerate via PanelPlot.export_html().
"""
import base64, io
import pandas as pd
import panel as pn
import param
import plotly.express as px

pn.extension("plotly")


#_data = pd.read_parquet(io.BytesIO(base64.b64decode(_DATA_B64.strip())))
_data = pd.read_parquet("{datapath}")

# ── PlotState ─────────────────────────────────────────────────────────────────
{state_src}

# ── Analysis ──────────────────────────────────────────────────────────────────
{analysis_src}

# ── Plot ──────────────────────────────────────────────────────────────────────
{plot_src}

# ── PanelPlot (minimal, self-contained, no base class needed) ─────────────────
class _App({state_cls_name}):
    def __init__(self, data, analysis, plot):
        super().__init__(
    {defaults_repr}
        )
        self._data     = data
        self._analysis = analysis
        self._plot     = plot

    def _param_names(self):
        return [n for n in {state_cls_name}.param if n != "name"]

    @property
    def _figure_pane(self):
        def _render(**kwargs):
            df  = self._analysis.run(self._data, self)
            fig = self._plot.render(df, self)
            return pn.pane.Plotly(fig, sizing_mode="stretch_width", height=520)
        return pn.bind(_render, **{{n: getattr(self.param, n) for n in self._param_names()}})

    def _sidebar(self):
        names  = self._param_names()
        groups = {{"Analysis": [], "Plot Layout": [], "Other": []}}
        for n in names:
            if n.startswith("analysis_"):   groups["Analysis"].append(n)
            elif n.startswith("plot_"):     groups["Plot Layout"].append(n)
            else:                           groups["Other"].append(n)
        sections = [(lbl, pn.Param(self, parameters=ns, show_name=False))
                    for lbl, ns in groups.items() if ns]
        return pn.Accordion(*sections, active=[0])

    def view(self):
        return pn.template.FastListTemplate(
            title="Interactive Analysis",
            sidebar=[self._sidebar()],
            main=[self._figure_pane],
        )

# ── entry point ───────────────────────────────────────────────────────────────
_app = _App(_data, {analysis_cls_name}(), {plot_cls_name}())
_app.view().servable()
'''


# Interactive Bioinformatics Plots → Standalone HTML

# Architecture

# PlotState          ← all tunable parameters (flat param.Parameterized)
#     ↓
# BaseAnalysis       ← pure computation: data + state → tidy DataFrame
#     ↓
# BasePlot           ← pure rendering:   tidy df + state → Plotly figure
#     ↓
# PanelPlot          ← UI + reactivity + WASM export

# PanelPlot.from_state() creates a dynamic subclass that inherits your
# PlotState directly, so all params are flat and pn.bind works without
# indirection.

# Naming convention for automatic widget grouping

# Prefix your param names:

# PrefixSidebar sectionanalysis_*Analysisplot_*Plot Layout(anything)Other

# Usage

# Interactive (local Panel server)

# pythonfrom example_pca import build_app, make_demo_data

# app = build_app(make_demo_data())
# app.show()           # opens browser at localhost:5006

# Export to standalone HTML

# pythonapp.export_html(outdir="dist/", stem="pca")
# # → dist/dist/pca.html   (share this file directly)

# Or from the command line:

# bashpython example_pca.py --export --outdir dist/

# What happens under the hood:


# Data is serialised to parquet and base64-encoded — embedded directly in app.py
# Your PlotState, Analysis, and Plot class sources are extracted via
# inspect.getsource()
# A self-contained app.py is rendered from a template (no mmalignments imports)
# panel convert app.py --to pyodide-worker --out dist/ is called
# Result: a single pca.html that runs entirely in the browser via WebAssembly


# Requirements

# panel>=1.3
# plotly
# pandas
# param
# pyarrow
# scikit-learn      # only for PCA example

# Adding a new plot type


# Subclass PlotState with your params
# Subclass BaseAnalysis and implement run(data, state) → DataFrame
# Subclass BasePlot and implement render(df, state) → plotly.Figure
# Call PanelPlot.from_state(MyState, data, MyAnalysis(), MyPlot())


# That's it. Export works automatically.

# Limitations of pyodide-worker export


# All data must be embedded (no network calls from the worker thread)
# Pure-Python packages only — packages with C extensions must be available in
# pyodide's package index (pandas, numpy, scikit-learn, plotly all are ✓)
# mmalignments itself cannot run in pyodide — the export embeds only the
# PlotState/Analysis/Plot classes, not your pipeline infrastructure
