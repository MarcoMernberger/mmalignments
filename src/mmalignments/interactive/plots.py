class Plot:

    def render(self, df, state: PlotState):
        import plotly.express as px

        return px.scatter(
            df,
            x="x",
            y="y",
            size=state.marker_size,
            color=state.face_color
        )