"""
app.py
───────
PlotConfig -> running Panel app. Your DAG node calls build_app() and gets
back something with .panel() (live) and .export_html() (standalone).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import panel as pn
import param
from param import Parameterized

from typing import Mapping

from .spec import PlotConfig, PlotLayer
from .plots import PLOT_REGISTRY, TitleState, ExportState
# fmt: off
from .orchestration import (
    build_combined_state_cls, 
    CompositeRenderer, 
    SelectiveRenderer, 
    build_selection_state_cls, 
    populate_dynamic_selectors, 
    build_color_maps
)
# fmt: on
pn.extension("plotly")


class InteractiveApp(Parameterized):
    """
    NOTE: this must itself be Parameterized — not a plain class — even
    though it has no params of its own. Two independent reasons:
      1. @param.depends / pn.bind's internal ref-resolution machinery
         calls `self.param` on whatever object owns the bound method;
         a plain class has no `.param` and raises AttributeError as soon
         as Panel tries to resolve dependencies (e.g. inside a Template).
      2. The actual dependency names (state_cls's params) are only known
         at RUNTIME, built dynamically per PlotConfig — @param.depends
         needs its argument names at decoration time, so we use pn.bind
         in the property below instead, binding directly to self.state's
         params once they exist.
    """

    DEFAULT_EXTRA_STATE = [TitleState]

    def __init__(
        self,
        data: pd.DataFrame | Mapping[str, pd.DataFrame],
        config: PlotConfig,
        title: str = "Interactive Analysis",
        **params,
    ):
        """ """
        self.title = title
        super().__init__(**params)
        if isinstance(data, pd.DataFrame):
            self.data = {"raw": data}
        else:
            self.data = data
        self.config = config
        self.select_param = config.select_param

        config.validate_against(PLOT_REGISTRY)
        build_color_maps(
            config.layers, self.data, explicit_maps=config.explicit_color_maps
        )

        extra_state = config.extra_state or self.DEFAULT_EXTRA_STATE

        # build the renderer
        if self.select_param is not None:
            select_state = build_selection_state_cls(
                config.layers, pname=self.select_param
            )
            extra_state.append(select_state)
            self.renderer = SelectiveRenderer(
                config.layers, select_param=self.select_param
            )
        else:
            self.renderer = CompositeRenderer(config.layers)

        # build combined class state
        self.state_cls = build_combined_state_cls(config.layers, extra_state)
        data_for_dynamic = self.data["raw"]
        populate_dynamic_selectors(
            self.state_cls,
            data_for_dynamic,
            result=None,
            shared_roles=config.shared_roles,
        )
        self.state = self.state_cls()

        self._last_fig = None
        self._last_result = None
        self._result_populated_from_analysis = False
        self._bind_actions()

    # ── reactive figure ───────────────────────────────────────────────────────

    @property
    def _figure_pane(self):
        """
        pn.bind, not @param.depends: the dependency set (every param in
        self.state_cls) is only known once state_cls is built at runtime,
        so it can't be written as a static decorator argument list.
        """
        param_names = [n for n in self.state_cls.param if n != "name"]
        return pn.bind(
            self._render, **{n: getattr(self.state.param, n) for n in param_names}
        )

    def _render(self, **kwargs):
        result = self.config.analysis.run(self.data, self.state)

        # Specs with from_result=True can only be resolved once we have
        # an AnalysisResult. Done once (not on every re-render) since
        # objects rarely change run-to-run; cheap guard either way.
        if not self._result_populated_from_analysis:
            populate_dynamic_selectors(
                self.state_cls,
                self.data["raw"],
                result=result,
                shared_roles=self.config.shared_roles,
            )
            self._result_populated_from_analysis = True

        self.config.validate_output(result, PLOT_REGISTRY)
        fig = self.renderer.render(result, self.state)
        self._last_fig = fig
        self._last_result = result
        return pn.pane.Plotly(fig, sizing_mode="stretch_width", height=520)

    # ── actions ───────────────────────────────────────────────────────────────

    def _bind_actions(self) -> None:
        self._handlers = {}
        for name in self.state_cls.param:
            if name == "name":
                continue
            if isinstance(self.state_cls.param[name], param.Action):
                self._handlers[name] = getattr(self, f"_action_{name}", self._noop)

    def _noop(self):
        pass

    def _action_export_png(self):
        self._export("png")

    def _action_export_svg(self):
        self._export("svg")

    def _action_export_tables(self):
        self._export_tables()

    def _action_export_app(self) -> None:
        """
        Generate a self-contained standalone_app.py (+ requirements.txt) that
        runs without the mmalignments package.

        Uses a lazy relative import so that:
          - in the full package, export.py is imported at call time (avoids a
            circular dependency at module load: export.py imports app.py).
          - in a standalone script (which doesn't ship the export module), the
            ImportError / SystemError is caught silently and the button is
            harmlessly inert.
        """
        try:
            from . import export as _exp
        except (ImportError, SystemError):
            print("[export_app] Export module not available in this environment.")
            return
        folder = getattr(self.state, "export_folder", "results")
        stem = getattr(self.state, "export_filename", "standalone_app")
        _exp.export_standalone_script(
            data=self.data,
            config=self.config,
            outdir=folder,
            stem=stem,
            title=self.title,
        )

    def _export(self, fmt: str) -> None:
        if self._last_fig is None:
            return
        folder = getattr(self.state, "export_folder", "results")
        filename = getattr(self.state, "export_filename", "figure")
        path = Path(folder) / f"{filename}.{fmt}"
        path.parent.mkdir(parents=True, exist_ok=True)
        self._last_fig.write_image(str(path))
        print(f"Saved {path}")

    def _export_tables(self) -> None:
        if self._last_result is None:
            return
        folder = getattr(self.state, "export_folder", "results")
        filename = getattr(self.state, "export_filename", "figure")
        path = Path(folder)
        path.parent.mkdir(parents=True, exist_ok=True)
        for suffix, df in self._last_result.data.items():
            df.to_csv(path / f"{filename}.{suffix}.tsv", index=False, sep="\t")
        print(f"Saved tables to {path}")

    # ── layout ───────────────────────────────────────────────────────────────

    def _action_param_names(self) -> set[str]:
        return {
            n
            for n in self.state_cls.param
            if isinstance(self.state_cls.param[n], param.Action)
        }

    def panel(self) -> pn.viewable.Viewable:
        action_names = self._action_param_names()
        all_names = [n for n in self.state_cls.param if n != "name"]

        # group by layer alias (each layer's namespaced params together),
        # everything else (TitleState/ExportState/the selector itself)
        # goes into "General"
        alias_prefixes = [layer.alias for layer in self.config.layers]
        groups: dict[str, list[str]] = {}
        for n in all_names:
            # matched = next((a for a in alias_prefixes if n.startswith(f"{a}_")), None)
            # group = matched or "General"
            group = self.state_cls.state_groups.get(n, "General")  # sidebar group
            groups.setdefault(group, []).append(n)

        sections = []
        group_order = self.config.sidebar_order or list(groups.keys())
        group_order = (
            group_order[: group_order.index("__LAYERS__")]
            + alias_prefixes
            + group_order[group_order.index("__LAYERS__") + 1 :]
        )
        for group_name in group_order:
            names = groups.get(group_name, [])
            normal = [n for n in names if n not in action_names]
            actions = [n for n in names if n in action_names]
            items = []
            if normal:
                items.append(pn.Param(self.state, parameters=normal, show_name=False))
            for a in actions:
                label = self.state_cls.param[a].label or a
                btn = pn.widgets.Button(name=label, button_type="success")
                btn.on_click(lambda e, n=a: self._handlers[n]())
                items.append(btn)

            section_content = pn.Column(*items)
            if self.select_param is not None and group_name in alias_prefixes:
                # Layer-specific section: only show its widgets while THIS
                # layer is the one currently selected. visible is bound
                # reactively to the selector param, so toggling the
                # dropdown shows/hides sections without rebuilding them.
                section_content = pn.Column(
                    *items,
                    visible=pn.bind(
                        lambda current, this_alias=group_name: current == this_alias,
                        getattr(self.state.param, self.select_param),
                    ),
                )

            sections.append((group_name, section_content))

        sidebar = pn.Accordion(*sections, active=[0])
        return pn.template.FastListTemplate(
            title=self.title,
            sidebar=[sidebar],
            main=[self._figure_pane],
            **self.config.panel_kwargs,
        )

    def show(self, **kwargs):
        self.panel().show(**kwargs)


def build_app(
    data: pd.DataFrame | Mapping[str, pd.DataFrame],
    config: PlotConfig,
    title: str = "Interactive Analysis",
) -> InteractiveApp:
    return InteractiveApp(data, config, title=title)
