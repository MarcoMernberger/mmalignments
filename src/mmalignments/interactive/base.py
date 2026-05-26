import param
import panel as pn
import io
from abc import abstractmethod

class BaseAnalysis(param.Parameterized):
    # Jede Unterklasse deklariert ihre eigenen Parameter,
    # Panel generiert daraus automatisch Widgets

    @abstractmethod
    def plot(self):
        """return ."""
        ...

    def panel(self):
        """Rendert Widgets + Plot nebeneinander."""
        widgets = pn.Param(self.param, show_name=False)
        plot_pane = pn.bind(lambda **kw: self.plot(), **{
            p: getattr(self.param, p) for p in self.param
            if p != 'name'
        })
        return pn.Row(widgets, pn.pane.Plotly(plot_pane))

    def export_png(self, path):
        self.plot().write_image(path, format="png")

    def export_svg(self, path):
        self.plot().write_image(path, format="svg")