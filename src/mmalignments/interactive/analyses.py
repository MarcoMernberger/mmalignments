class Analysis:

    def run(self, data, state: PlotState):
        # ONLY compute, no UI
        df = data.copy()

        if state.normalize:
            df["value"] = df["value"] / df["value"].max()

        return df