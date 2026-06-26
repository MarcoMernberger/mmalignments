from .spec import Analysis, AnalysisResult
from pandas import DataFrame
from typing import Mapping


class PassthroughAnalysis(Analysis):
    def run(self, data: Mapping[str, DataFrame], state) -> AnalysisResult:
        df = data["raw"].copy()
        if getattr(state, "analysis_treatment", None):
            df = df[df["treatment"] == state.analysis_treatment]
        if getattr(state, "analysis_filter", ""):
            df = df.query(state.analysis_filter)
        return AnalysisResult(raw=df, agg=None)
