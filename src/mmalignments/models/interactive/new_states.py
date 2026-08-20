import param

from .new_specs import PlotState, PlotStateSpec

# registry for plt features

FEATURE_PARAMS = {
    "opacity": {
        "opacity": param.Number(default=1, bounds=(0, 1)),
        "line_width": param.Number(default=1),
    },
    "axelabels": {
        "x_label": param.String(""),
        "y_label": param.String(""),
    },
    "title": {
        "title": param.String(""),
    },
    "legend": {
        "show_legend": param.Boolean(True),
    },
}


class PlotStateFactory:
    """A factory to dynamically create PlotState objects based on a plot specification. This allows for flexible and reusable plot configurations."""

    def build(self, spec: PlotStateSpec):

        params = {}

        # shared features
        for feature in spec.features:
            params.update(FEATURE_PARAMS[feature])

        # plot-specific
        params.update(spec.params)

        # dynamic role selectors
        for name, role in spec.roles.items():
            if role.source == "parameter":
                params[name] = param.Selector(
                    default=role.columns[0],
                    objects=list(role.columns),
                )
            elif role.source == "dynamic":
                params[name] = param.Selector(
                    default=None,
                    objects=[None],
                    label=role,
                )
            elif role.source == "column":
                params[name] = param.Selector(
                    default=role.columns[0],
                    objects=list(role.columns),
                    label=role,
                )
        params["sidebar_group"] = spec.alias

        return type(f"{spec.alias}PlotState", (PlotState,), params)


# usage HBarPlotState = PlotStateFactory().build(HBarPlotSpec)
