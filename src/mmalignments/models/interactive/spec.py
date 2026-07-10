"""
spec.py
────────
Core contracts. Fixes vs. the draft:
  - PlotLayer is data, plot classes are registered separately (registry.py)
  - roles is consistently a dict[str, str] everywhere (role -> column name)
  - AnalysisResult keeps raw/agg as named fields (type-checkable), with both
    dict-style __getitem__ AND attribute access so layer.source resolves
    cleanly either way
  - Analysis ABC signature is consistent: run(data, state) -> AnalysisResult
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Mapping
from unicodedata import name

import plotly.graph_objects as go
import pandas as pd  # type: ignore[import]
import param  # type: ignore[import]
from param import Parameterized  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

# ─────────────────────────────────────────────────────────────────────────────
# AnalysisResult — fixed fields (raw/agg), but subscriptable by name too
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class AnalysisResult:
    """
    data: Mapping[str, pd.DataFrame] = {"raw": raw_df, "agg": agg_df}
    raw : always the per-datapoint frame (one row per observation)
    agg : aggregated frame, or None if this analysis never aggregates

    Subscript access (result["raw"], result["agg"]) lets PlotLayer.source
    stay a plain string without the renderer needing to know about
    dataclass field access vs dict access.
    """

    data: Mapping[str, pd.DataFrame] = field(default_factory=dict)

    def __getitem__(self, key: str) -> pd.DataFrame:
        if key not in self.data:
            raise KeyError(
                f"AnalysisResult has no source '{key}'. Available: {list(self.data.keys())}"
            )
        value = self.data[key]
        if value is None:
            raise KeyError(
                f"AnalysisResult.{key} is None — this analysis didn't "
                f"produce an aggregated frame. Check your PlotLayer.source."
            )
        return value


class Analysis(ABC):
    """data + state -> AnalysisResult. state is whatever Parameterized
    instance InteractiveApp builds (the combined state across all layers)."""

    @abstractmethod
    def run(
        self, data: Mapping[str, DataFrame], state: Parameterized
    ) -> AnalysisResult: ...


# ─────────────────────────────────────────────────────────────────────────────
# RoleSpec — the DATA contract a plot class declares (separate from UI knobs)
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class RoleSpec:
    name: str  # "x", "y", "color"
    dtype: str = "any"  # "continuous", "categorical", "any"
    required: bool = True
    description: str = ""


# ─────────────────────────────────────────────────────────────────────────────
# PlotLayer — ONE instance of a plot type in the composite figure.
# Namespacing key for Problem 1 is `alias`, not `plot_type`.
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class PlotLayer:
    plot_type: str  # registry key, e.g. "rate", "hbar"
    alias: str  # unique id for THIS instance
    roles: dict[str, str]  # role name -> actual column name
    source: str = "raw"  # which AnalysisResult field to read
    hover: list[str] = field(default_factory=list)
    color_maps: dict[str, dict] = field(default_factory=dict)


# ─────────────────────────────────────────────────────────────────────────────
# PlotConfig — single declarative source of truth for one app
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class PlotConfig:
    analysis: Analysis  # the Analysis class to use
    layers: list[PlotLayer]  # Plot Layers to render
    labels: dict[str, str] = field(default_factory=dict)
    extra_state: list[type] = field(default_factory=list)
    extra_state_bottom: list[type] = field(default_factory=list)
    select_param: str | None = None  # optional, only needed for SelectiveRenderer
    shared_roles: dict[str, str] = field(
        default_factory=dict
    )  # role -> column name, used by DynamicSelectorSpec.role resolution
    sidebar_order: list[str] | None = None  # optional, order of sidebar groups
    panel_kwargs: dict[str, Any] = field(
        default_factory=dict
    )  # additional kwargs for the panel template
    explicit_color_maps: dict[str, dict[str, str]] = field(
        default_factory=dict
    )  # override dynamic coloring for specific layers (e.g. {"count": {"Rez": "#1f77b4", "DMSO": "#ff7f0e"}})

    def validate_against(self, registry: dict[str, type]) -> None:
        """Check every layer's roles satisfy its plot class's REQUIRED roles."""
        problems = []
        for layer in self.layers:
            if layer.plot_type not in registry:
                problems.append(
                    f"'{layer.alias}': unknown plot_type '{layer.plot_type}'"
                )
                continue
            plot_cls = registry[layer.plot_type]
            missing = [r for r in plot_cls.required_roles() if r not in layer.roles]
            if missing:
                problems.append(
                    f"'{layer.alias}' ({layer.plot_type}) missing roles: {missing}"
                )
        if problems:
            raise ValueError("PlotConfig invalid:\n  " + "\n  ".join(problems))

    def validate_output(
        self,
        result: AnalysisResult,
        registry: dict[str, type],
        state: Parameterized | None = None,
    ) -> None:
        """After Analysis.run(): check the promised columns actually exist."""
        for layer in self.layers:
            df = result[layer.source]
            plot_cls = registry[layer.plot_type]
            plot_cls.validate_roles(df, layer.roles, state)


# ─────────────────────────────────────────────────────────────────────────────
# PlotState — base for per-plot-type LAYOUT knobs (UI contract, not data)
# ─────────────────────────────────────────────────────────────────────────────


class PlotState(Parameterized):
    """Subclass per plot type for its layout parameters only
    (marker_size, linewidth, ...). Never put role/column mappings here —
    those live in PlotLayer.roles.

    sidebar_group : where this state's widgets appear in the app sidebar.
        Defaults to the class name if not overridden. Layer-specific
        states (state_cls of a BasePlot) are still namespaced/grouped by
        layer ALIAS regardless of this — this only matters for extra_state
        classes (TitleState, ExportState, your own cross-cutting states),
        which otherwise all fall into one undifferentiated "General" bucket.
    """

    sidebar_group: str = "General"


# ─────────────────────────────────────────────────────────────────────────────
# DynamicSelectorSpec — declares a role-driven, data-populated selector
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class DynamicSelectorSpec:
    """
    Declares that one param on a PlotState subclass must have its
    `.objects` populated from real data at app-build time, rather than
    being hardcoded — analogous to how PlotLayer.roles maps a generic
    role name to an actual column, except here the role drives a
    Selector's choices instead of an axis.

    param_name : the param on the state class to populate (e.g. "analysis_treatment")
    role       : generic role name, resolved against `roles` at build time
                 (e.g. "treatment" -> "treatment_col" if your data uses a
                 differently-named column — same indirection PlotLayer uses)
    source     : "raw" to read column values straight from the input
                 DataFrame (before Analysis.run()), or "agg"/"raw" to read
                 from the AnalysisResult instead (after Analysis.run()) —
                 set via `from_result=True` when the column only exists
                 post-analysis (e.g. a computed "enrichment_bin" column)
    from_result: if True, source is resolved against AnalysisResult[source]
                 instead of the raw input DataFrame.
    include_none: prepend None to the resulting objects list (for an
                 "all / no filter" choice) — matches your
                 TreatmentSelectionState's objects=[None, ...] pattern.
    """

    param_name: str
    role: str
    source: str = "raw"
    from_result: bool = False
    include_none: bool = True


class DynamicPlotState(PlotState):
    """
    Base for PlotState subclasses that need data-driven Selector objects.
    Subclasses declare DYNAMIC_SELECTORS; the app-builder calls
    populate_dynamic_selectors() once data/roles are known, before the
    state instance is exposed to the UI.

        class TreatmentSelectionState(DynamicPlotState):
            DYNAMIC_SELECTORS = [
                DynamicSelectorSpec(param_name="analysis_treatment", role="treatment"),
            ]
            analysis_treatment = param.Selector(default=None, objects=[None], label="Treatment")
    """

    DYNAMIC_SELECTORS: list[DynamicSelectorSpec] = []


# ─────────────────────────────────────────────────────────────────────────────
# BasePlot — the plot class contract
# ─────────────────────────────────────────────────────────────────────────────


class BasePlot:
    ROLE_SPECS: dict[str, RoleSpec] = {}
    state_cls: type[PlotState] = PlotState
    trace_capabilities: set[str] = set()

    @classmethod
    def get_param_from_roles(cls, name: str, state: PlotState | None = None) -> str:
        if state is not None:
            return getattr(state, name, name)
        else:
            return getattr(cls.state_cls, name, name)

    @classmethod
    def required_roles(cls) -> set[str]:
        return {name for name, spec in cls.ROLE_SPECS.items() if spec.required}

    @classmethod
    def optional_roles(cls) -> set[str]:
        return {name for name, spec in cls.ROLE_SPECS.items() if not spec.required}

    @classmethod
    def validate_roles(
        cls, df: pd.DataFrame, roles: dict[str, str], state: PlotState | None = None
    ) -> None:
        for role_name, spec in cls.ROLE_SPECS.items():
            if role_name not in roles:
                if spec.required:
                    raise ValueError(
                        f"{cls.__name__}: missing required role '{role_name}'"
                    )
                continue
            col = cls.get_param_from_roles(roles[role_name], state=state)
            if col not in df.columns:
                raise ValueError(
                    f"{cls.__name__}: role '{role_name}' maps to column "
                    f"'{col}', which is not in the DataFrame "
                    f"(columns: {list(df.columns)}) or the PlotState: {dir(state)}"
                )

    def render(
        self,
        df: pd.DataFrame,
        layer: PlotLayer,
        state: PlotState,
        processor: "TraceProcessor",
        fig=None,
    ) -> go.Figure:
        raise NotImplementedError
