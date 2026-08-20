from __future__ import annotations

from dataclasses import dataclass, import field
from typing import Literal


@dataclass
class RoleSpec:
    name: str
    dtype: str
    required: bool = True

    source: Literal["column", "parameter", "dynamic"] = "column"

    columns: tuple[str, ...] = ()

    parameter: str | None = None

    dynamic_column: str | None = None


@dataclass
class PlotLayer:
    plot_type: str  # registry key, e.g. "rate", "hbar"
    alias: str  # unique id for THIS instance
    roles: dict[str, RoleSpec]  # role name -> RoleSpec
    source: str = "raw"  # which AnalysisResult field to read
    hover: list[str] = field(default_factory=list)
    color_maps: dict[str, dict] = field(default_factory=dict)


# example: fixed dataframe colum
# RoleSpec(
#     name="x",
#     dtype="continuous",
#     source="column",
#     columns=("count",)
# )
# this would mean, that x is to be taken from the column "count" in the dataframe, and it is of type continuous.

# example: selector for multiple columns
# RoleSpec(
#     name="x",
#     dtype="continuous",
#     source="parameter",
#     parameter="score",
#     columns=(
#         "count",
#         "logFC",
#         "beta",
#     )
# )
# this would mean, that x is to be taken from the column specified by the parameter "score", which can be one of "count", "logFC", or "beta", and it is of type continuous.

# example: dynamic selector from column uniqeue values
# RoleSpec(
#     name="group",
#     dtype="categorical",
#     source="dynamic",
#     dynamic_column="condition"
# )


@dataclass
class PlotStateSpec:
    alias: str
    roles: dict[str, RoleSpec]
    features: list[str] = field(default_factory=list)
    # state_cls: type[PlotState]
