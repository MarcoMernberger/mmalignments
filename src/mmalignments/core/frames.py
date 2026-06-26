from typing import Callable

from pandas import (  # a couple of reusable dataframe functions often used
    DataFrame,
    Series,
)


def append(df: DataFrame, call: Callable[[], tuple[str, Series]],) -> DataFrame:
    """Add a column to the dataframe based on the provided scoring function."""
    name, series = call()
    df[name] = series
    return df

def drop(df: DataFrame, *columns: str) -> DataFrame:
    """Drop specified columns from the dataframe."""
    df = df.drop(columns=columns)
    return df

def ensure_types(df: DataFrame, scheme: dict[str, type]) -> DataFrame:
    """Ensure the dataframe columns have the specified types."""
    for column, dtype in scheme.items():
        df[column] = df[column].astype(dtype)
    return df

def merge_frames(df: DataFrame, other: DataFrame, **kwargs) -> DataFrame:
    """Merge two dataframes on specified keys."""
    left_on = kwargs.get("left_on", "Flanks")
    right_on = kwargs.get("right_on", "Annotation")
    how = kwargs.get("how", "inner")
    df = df.merge(other, left_on=left_on, right_on=right_on, how=how)
    return df


def group(df: DataFrame, key_columms: list[str], aggregation: dict[str, str]) -> DataFrame:
    df = (
        df.groupby(
            key_columms,
            as_index=False,
        )
        .agg(
            aggregation
        )
    )
    return df


def groupmelt(df: DataFrame, id_vars: list[str], value_vars: list[str], var_name: str, value_name: str) -> DataFrame:
    """Convert a dataframe from wide to long format."""
    df.index = pd.MultiIndex.from_frame(df[id_vars])

    df = df.melt(
        id_vars=id_vars,
        value_vars=value_vars,
        var_name=var_name,
        value_name=value_name
    )
    return df
