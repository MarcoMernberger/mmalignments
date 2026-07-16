"""
orchestration.py
──────────────────
Combines N PlotLayers into one Parameterized state class, namespaced by
ALIAS (not plot_type) — this is what solves Problem 1: two layers of the
same plot_type (or different ones) never collide, because each layer
already carries a unique alias by construction.
"""

from __future__ import annotations

from typing import Mapping

import pandas as pd  # type: ignore[import]
import param  # type: ignore[import]
from param import Parameterized  # type: ignore[import]

from .plots import PLOT_REGISTRY, TraceProcessor
from .spec import (
    AnalysisResult,
    DynamicPlotState,
    DynamicSelectorSpec,
    LayerStateView,
    PlotLayer,
    PlotState,
)

# Same ordering Plotly Express uses by default, so "auto" mode looks
# identical to what you'd get without any color_map — just stable across
# filtering instead of recomputed from whatever's left after a filter.
_DEFAULT_PALETTE = [
    "#636EFA",
    "#EF553B",
    "#00CC96",
    "#AB63FA",
    "#FFA15A",
    "#19D3F3",
    "#FF6692",
    "#B6E880",
    "#FF97FF",
    "#FECB52",
]

################################################################################
# Custom Color Maps
################################################################################


def build_color_maps(
    layers: list[PlotLayer],
    data: Mapping[str, pd.DataFrame],
    explicit_maps: dict[str, dict[str, dict]] | None = None,
) -> None:
    """
    Populates layer.color_maps["color"] for every layer whose roles include
    a "color" role, using the FULL (unfiltered) category set from `data` —
    so Plotly always assigns the same color to the same category, no
    matter how Analysis.run() later filters the data per-render.

    Mutates each PlotLayer in place (color_maps is a plain dict field, not
    a param — safe to mutate directly, unlike the Selector .objects case).

    explicit_maps : {layer_alias: {category_value: color}} — for layers
        where you want a SPECIFIC fixed mapping (e.g. {"Rez": "#1f77b4",
        "DMSO": "#ff7f0e"}) instead of the auto-generated palette. Layers
        not mentioned here fall back to automatic assignment from the
        full category set, in sorted order (deterministic, not
        insertion-order-dependent — insertion order would itself be
        fragile across different filtered subsets).

    Call this ONCE, right after loading the input DataFrame and before
    building the combined state / first render — analogous to
    populate_dynamic_selectors() for Selector objects, except color_maps
    live on PlotLayer (read every render by the plot's render()), not on
    a param.
    """
    explicit_maps = explicit_maps or {}
    for layer in layers:
        if "color" not in layer.roles:
            continue
        column = layer.roles["color"]
        if column not in data[layer.source].columns:
            raise ValueError(
                f"PlotLayer '{layer.alias}': color role maps to column "
                f"'{column}', not found in data. Available: {list(data[layer.source].columns)}"
            )

        if layer.alias in explicit_maps:
            layer.color_maps["color"] = dict(explicit_maps[layer.alias])
            continue

        categories = sorted(
            data[layer.source][column].dropna().unique().tolist(), key=str
        )
        palette = _DEFAULT_PALETTE
        color_map = {cat: palette[i % len(palette)] for i, cat in enumerate(categories)}
        layer.color_maps["color"] = color_map


# ─────────────────────────────────────────────────────────────────────────────
# Parameter copying helper (param.Parameter instances are NOT shareable
# across classes — must be freed of their previous "owner" binding)
# ─────────────────────────────────────────────────────────────────────────────


def _fresh_param_copy(p: param.Parameter) -> param.Parameter:
    """
    param.Parameter instances are bound for life to their defining class
    (via the `name`/`owner` slots) and refuse to be re-bound — even
    copy.copy() carries that binding over and raises on first use in a
    new class. There is no public "unbind" API, so we reconstruct a fresh
    instance of the same type from its own slots instead of copying the
    bound object.

    We use _all_slots_ (includes base Parameter slots PLUS subclass-specific
    ones like `bounds`/`step` for Number, `_objects` for Selector, etc.) so
    this works generically for any param type without a type-specific
    branch per parameter class.

    CAUTION: do not blanket-skip every None value here. Several slots use
    None as a legitimate, meaningful value (e.g. a String param declared
    with default=None, allow_None=True — common for "unset" labels like
    appearance_xlabel). Skipping None unconditionally silently replaces
    such defaults with the Parameter subclass's OWN fallback default
    (e.g. "" for String), which is a different value, not a copy.

    The one slot that genuinely cannot be passed as explicit None in some
    param versions is `default_factory` — passing None there raises
    ("must be a callable, not NoneType") even though *omitting* it
    entirely is fine and equivalent. So that single slot is special-cased;
    every other slot — None or not — is forwarded as-is.

    `names` (on Selector) is excluded outright: it's a derived, read-only
    mapping (label -> value, populated automatically when `objects` is a
    dict) rather than a constructor argument — passing it raises
    "unexpected keyword argument 'names'".
    """
    cls = type(p)
    if isinstance(p, param.Action):
        return param.Action(
            default=p.default,
            label=p.label,
            doc=p.doc,
            precedence=p.precedence,
        )
    exclude = {"name", "owner", "watchers", "names", "length"}
    rename = {"_label": "label", "_objects": "objects"}
    kwargs: dict = {}
    for slot in cls._all_slots_:
        if slot in exclude:
            continue
        key = rename.get(slot, slot)
        try:
            value = getattr(p, key)
        except AttributeError:
            continue
        if key == "default_factory" and value is None:
            continue  # the one slot where None must be omitted, not passed
        kwargs[key] = value

    if isinstance(p, param.Selector):
        kwargs["objects"] = (
            p._objects
        )  # use the raw _objects, not the derived .objects property
    try:
        res = cls(**kwargs)
    except Exception as e:
        print(p.name)
        print(p)
        print(kwargs)
        raise e
    return res


################################################################################
# State Builder
################################################################################

# DynamicSelector States


def populate_dynamic_selectors(
    combined_cls: type[Parameterized],
    data: Mapping[str, pd.DataFrame],
    result: AnalysisResult | None,
    shared_roles: dict[str, str],
) -> None:
    """
    Resolves every DynamicSelectorSpec collected during
    build_combined_state_cls() against real data, setting `.objects` on
    the relevant Selector params IN PLACE on the class (so every instance
    of combined_cls picks them up).

    Must be called AFTER build_combined_state_cls() and — for specs with
    from_result=True — AFTER Analysis.run() has produced `result`. Specs
    with from_result=False only need `data` and can be resolved earlier
    (e.g. right after loading the input DataFrame, before the first render).

    Why mutate the class instead of the instance: param.Selector.objects
    must be set before any instance reads/validates against it, and the
    sidebar widget (pn.Param) is built from the CLASS's parameter
    definitions. Setting it on the class once, here, is the simplest
    correct point — equivalent to what your original draft's
    `init_from_frame` did per-instance, but the Selector's `objects` slot
    is inherently class-level metadata in param, not per-instance state.
    """
    for pname, dyn_spec in combined_cls._dynamic_specs:
        if dyn_spec.from_result and result is None:
            # Not an error: from_result specs are resolved later, once
            # Analysis.run() has produced a result (see InteractiveApp._render,
            # which calls this function again with result populated).
            # Skip silently on this pass instead of every spec needing the
            # caller to pre-filter by from_result.
            continue
        if dyn_spec.role not in shared_roles:
            raise ValueError(
                f"DynamicSelectorSpec for '{pname}' needs role "
                f"'{dyn_spec.role}', but PlotConfig.shared_roles has no "
                f"mapping for it. Add e.g. "
                f"shared_roles={{'{dyn_spec.role}': 'your_column_name'}}."
            )
        column = shared_roles[dyn_spec.role]
        if dyn_spec.from_result:
            if result is None:
                raise ValueError(
                    f"DynamicSelectorSpec for '{pname}' has from_result=True "
                    f"but no AnalysisResult was provided yet."
                )
            source_df = result[dyn_spec.source]
        else:
            source_df = data[dyn_spec.source]

        if column not in source_df.columns:
            raise ValueError(
                f"DynamicSelectorSpec for '{pname}': column '{column}' "
                f"(role '{dyn_spec.role}') not found in "
                f"{'AnalysisResult.' + dyn_spec.source if dyn_spec.from_result else 'input data'}. "
                f"Available columns: {list(source_df.columns)}"
            )
        values = sorted(source_df[column].dropna().unique().tolist())
        objects = dyn_spec.objects
        names = dyn_spec.names
        if isinstance(dyn_spec.objects, list):
            objects = dyn_spec.objects + values
            if dyn_spec.names:
                names = dyn_spec.names + [str(v) for v in values]
                objects = dict(zip(names, objects))  # convert to dict for param.Selector
        else:
            objects = ([None] + values) if dyn_spec.include_none else values
        param_obj = combined_cls.param[pname]
        if not isinstance(param_obj, param.Selector):
            raise TypeError(
                f"DynamicSelectorSpec targets '{pname}', but that param is "
                f"a {type(param_obj).__name__}, not a Selector."
            )
        param_obj.objects = objects # this normalizes objects to list and sets "names"
        # keep default valid against the new objects list
        if param_obj.default not in objects:
            if isinstance(param_obj, param.ListSelector):
                if isinstance(objects, dict):
                    param_obj.default = objects.values()
                else:
                    param_obj.default = objects
            elif isinstance(param_obj, param.Selector):
                if isinstance(objects, dict):
                    param_obj.default = next(iter(objects.values()))
                else:
                    param_obj.default = objects[0]
            else:
                pass

# Combined state builder — namespaced by LAYER ALIAS


def build_combined_state_cls(
    layers: list[PlotLayer],
    extra_state: list[type[PlotState]] = (),
) -> type[Parameterized]:
    """
    For layers=[PlotLayer(alias="count", plot_type="rate"),
                PlotLayer(alias="frequency", plot_type="rate")]
    produces param names:
        count_marker_size, count_linewidth, ...
        frequency_marker_size, frequency_linewidth, ...
    Even though both layers use the SAME plot_type, they never collide
    because namespacing is by alias, not plot_type.

    extra_state params are merged in unprefixed, shared across all layers.
    Also records, per resulting class, which extra_state CLASS each
    unprefixed param came from (_extra_state_source) and that class's
    sidebar_group — this is what lets the app's panel() put each
    extra_state param into its own named sidebar section instead of one
    undifferentiated "General" bucket.
    """
    namespace: dict[str, param.Parameter] = {}
    layer_param_map: dict[str, dict[str, str]] = {}
    # pname -> (sidebar_group, source class) for extra_state params only
    state_groups: dict[str, str] = {}
    dynamic_specs: list[tuple[str, DynamicSelectorSpec]] = []  # (pname, spec)
    for layer in layers:
        plot_cls = PLOT_REGISTRY[layer.plot_type]
        state_cls = plot_cls.state_cls
        own_map: dict[str, str] = {}
        for pname in state_cls.param:
            if pname == "name":
                continue
            ns_name = f"{layer.alias}_{pname}"
            namespace[ns_name] = _fresh_param_copy(state_cls.param[pname])
            own_map[pname] = ns_name
            state_groups[ns_name] = layer.alias
        layer_param_map[layer.alias] = own_map
        if issubclass(state_cls, DynamicPlotState):
            for dyn_spec in state_cls.DYNAMIC_SELECTORS:
                dynamic_specs.append((f"{layer.alias}_{dyn_spec.param_name}", dyn_spec))

    for extra_cls in extra_state:
        group = getattr(extra_cls, "sidebar_group", "General")
        own_map: dict[str, str] = {}
        for pname in extra_cls.param:
            if pname == "name" or pname in namespace:
                continue
            ns_name = f"{group}_{pname}"
            namespace[ns_name] = _fresh_param_copy(extra_cls.param[pname])
            # print(f"Adding extra_state param {pname} as {ns_name} in group {group}")
            state_groups[ns_name] = group
            own_map[pname] = ns_name

        if issubclass(extra_cls, DynamicPlotState):
            for dyn_spec in extra_cls.DYNAMIC_SELECTORS:
                dynamic_specs.append((f"{group}_{dyn_spec.param_name}", dyn_spec))
        if group in layer_param_map:
            layer_param_map[group].update(own_map)
        else:
            layer_param_map[group] = own_map

    combined_cls = type("CombinedState", (Parameterized,), namespace)
    combined_cls._layer_param_map = layer_param_map
    combined_cls.state_groups = state_groups
    combined_cls._dynamic_specs = dynamic_specs
    return combined_cls


# SelectionState builder — a tiny PlotState subclass containing ONLY


def build_layer_selector_param(
    layers: list[PlotLayer],
    pname: str = "analysis_layer_select",
    label: str = "View",
) -> param.Selector:
    """
    Builds a Selector whose `objects` are the ACTUAL layer aliases from
    this config — not a hardcoded list. This replaces a static
    SelectionState class (which can't know your aliases in advance) with
    a factory that derives the selector from config.layers at the point
    build_combined_state_cls() is called.

    Use together with build_combined_state_cls() by passing a tiny
    one-off extra_state class built around this param — see
    build_selection_state_cls() below for the convenience wrapper.
    """
    aliases = [layer.alias for layer in layers]
    if not aliases:
        raise ValueError("build_layer_selector_param(): layers is empty")
    res = param.Selector(default=aliases[0], objects=aliases, label=label)
    return res


def build_selection_state_cls(
    layers: list[PlotLayer],
    pname: str = "analysis_layer_select",
    label: str = "View",
) -> type[PlotState]:
    """
    Convenience wrapper: builds a tiny PlotState subclass containing ONLY
    the layer-selector param, with `objects` derived from the real
    layers — pass this in extra_state instead of a hand-written
    SelectionState with hardcoded objects=[...].

        layers = [PlotLayer(plot_type="rate", alias="count", ...),
                   PlotLayer(plot_type="rate", alias="frequency", ...)]
        SelectionState = build_selection_state_cls(layers)
        state_cls = build_combined_state_cls(layers, extra_state=[SelectionState, TitleState])
    """
    selector = build_layer_selector_param(layers, pname=pname, label=label)
    return type(f"{pname}", (PlotState,), {pname: selector})


def build_dynam_state_cls(
    layers: list[PlotLayer],
    pname: str = "analysis_layer_select",
    label: str = "View",
) -> type[PlotState]:
    """
    Convenience wrapper: builds a tiny PlotState subclass containing ONLY
    the layer-selector param, with `objects` derived from the real
    layers — pass this in extra_state instead of a hand-written
    SelectionState with hardcoded objects=[...].

        layers = [PlotLayer(plot_type="rate", alias="count", ...),
                   PlotLayer(plot_type="rate", alias="frequency", ...)]
        SelectionState = build_selection_state_cls(layers)
        state_cls = build_combined_state_cls(layers, extra_state=[SelectionState, TitleState])
    """
    selector = build_layer_selector_param(layers, pname=pname, label=label)
    return type(f"{pname}", (PlotState,), {pname: selector})

#  deprecated ... no more double copying of state, just a view on the combined state
# def extract_layer_state(combined: Parameterized, layer: PlotLayer, plot_cls: type):
#     """Build a real, lightweight state_cls instance for ONE layer, populated
#     from the combined state's namespaced values. The plot class's render()
#     never sees aliases or namespacing — just its own plain state object."""
#     own_map = combined._layer_param_map[layer.alias]

#     print("Extracting layer state for", layer.alias)
#     print(combined._layer_param_map[layer.alias])
#     sub_state = plot_cls.state_cls()
#     for raw_name, ns_name in own_map.items():
#         setattr(sub_state, raw_name, getattr(combined, ns_name))
#     return sub_state


def construct_layer_state_view(combined: Parameterized, layer: PlotLayer):
    """Build a view on the combined state for a specific layer."""
    own_map = combined._layer_param_map[layer.alias]
    return LayerStateView(combined, own_map)

# ─────────────────────────────────────────────────────────────────────────────
# CompositeRenderer — dispatches EVERY layer to its plot class, merges figures
# ─────────────────────────────────────────────────────────────────────────────


class CompositeRenderer:
    def __init__(self, layers: list[PlotLayer]):
        self.layers = layers
        self.processor = TraceProcessor()

    def render(
        self, result: AnalysisResult, combined_state: Parameterized, layer_state_views: dict[str, LayerStateView]
    ) -> "go.Figure":
        fig = None
        for layer in self.layers:
            plot_cls = PLOT_REGISTRY[layer.plot_type]
            plot = plot_cls()
            df = result[layer.source]
            layer_state_view = layer_state_views[layer.alias]

            fig = plot.render(df, layer, layer_state_view, combined_state, self.processor, fig=fig)

        # Shared appearance (title, axis labels, tick relabeling) is applied
        # ONCE at the end, on the combined figure, reading directly from
        # combined_state — no grafting onto per-layer state needed.
        if fig is not None:
            fig = self.processor.apply_shared(fig, combined_state, self.layers, result)
        return fig


# ─────────────────────────────────────────────────────────────────────────────
# SelectiveRenderer — renders ONLY the layer currently chosen via a
# Selector param (built by build_layer_selector_param above), instead of
# combining all layers into one figure. Useful for "switch between views"
# UIs (e.g. toggle between Count / Frequency / Enrichment).
# ─────────────────────────────────────────────────────────────────────────────


class SelectiveRenderer:
    def __init__(
        self, layers: list[PlotLayer], select_param: str = "analysis_layer_select"
    ):
        """
        layers       : the full list of candidate layers (same list you'd
                        pass to build_combined_state_cls / build_selection_state_cls)
        select_param : the param name holding the user's current choice —
                        must match the pname used in build_selection_state_cls()
        """
        self.layers_by_alias: dict[str, PlotLayer] = {la.alias: la for la in layers}
        self.layers = layers  # kept for apply_shared's tick-label lookup
        self.select_param = select_param
        self.processor = TraceProcessor()

    def render(
        self, result: AnalysisResult, combined_state: Parameterized, layer_state_views: dict[str, LayerStateView]
    ) -> "go.Figure":
        alias = getattr(combined_state, self.select_param, None)
        if alias is None:
            raise ValueError(
                f"SelectiveRenderer.render(): combined_state has no "
                f"'{self.select_param}' param — did you forget to include "
                f"build_selection_state_cls(layers) in extra_state when "
                f"building the app?"
            )
        if alias not in self.layers_by_alias:
            raise ValueError(
                f"SelectiveRenderer.render(): '{alias}' is not a known layer "
                f"alias. Available: {list(self.layers_by_alias)}"
            )

        layer = self.layers_by_alias[alias]
        plot_cls = PLOT_REGISTRY[layer.plot_type]
        plot = plot_cls()
        df = result[layer.source]
        #layer_state = extract_layer_state(combined_state, layer, plot_cls)
        layer_state_view = layer_state_views[layer.alias]
        fig = plot.render(df, layer, layer_state_view, combined_state, self.processor, fig=None)
        if fig is not None:
            fig = self.processor.apply_shared(fig, combined_state, [layer], result)
        return fig
