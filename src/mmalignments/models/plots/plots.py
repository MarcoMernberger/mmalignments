from abc import ABC, abstractmethod

import panel as pn


class PanelPlot(ABC):

    def __init__(self):
        # initalize controls based on supported_controls
        self.controls = self._build_controls()

    @property
    def supported_controls(self) -> set[str]:
        """
        which adjustable plot params are supported by this plot? This determines which 
        controls are shown in the UI.

        _extended_summary_

        Returns
        -------
        set[str]
            _description_
        """
        return set()

    # -------------------------
    # 2. Standard Controls
    # -------------------------
    def _build_controls(self):
        controls = {}

        controls["line_width"] = pn.widgets.FloatSlider(
            name="Line width", start=0.5, end=5.0, value=1.5
        )

        controls["marker_size"] = pn.widgets.IntSlider(
            name="Marker size", start=5, end=200, value=40
        )

        controls["face_color"] = pn.widgets.Select(
            name="Face color",
            options=["steelblue", "tomato", "seagreen", "gray"],
            value="steelblue",
        )

        controls["line_color"] = pn.widgets.Select(
            name="Line color", options=["black", "red", "blue", "green"], value="black"
        )

        controls["font_size"] = pn.widgets.IntSlider(
            name="Font size", start=6, end=20, value=10
        )

        controls["xlabel"] = pn.widgets.TextInput(name="X label", value="X")
        controls["ylabel"] = pn.widgets.TextInput(name="Y label", value="Y")

        return {k: v for k, v in controls.items() if k in self.supported_controls}

    # -------------------------
    # 3. UI Assembly
    # -------------------------
    def view(self):
        bound_plot = pn.bind(self._plot, **self.controls)
        return pn.Column(pn.Row(*self.controls.values()), bound_plot)

    # -------------------------
    # 4. MUSS implementiert werden
    # -------------------------
    @abstractmethod
    def _plot(self, **kwargs):
        pass
