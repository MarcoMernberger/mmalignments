"Some helper functions for setting up a run of the mmalignments pipeline."

from pandas import DataFrame, Series  # type: ignore[import]


def samples_by_factors(samples: DataFrame, **factors) -> list[str]:
    """
    returns a dictionary of group names to sample names based on the provided
    factors.

    Parameters
    ----------
    samples : DataFrame
        The samples DataFrame containing sample and factor information.

    Returns
    -------
    list[str]
        A list of sample names that match the provided factors.
    """
    mask = Series([True] * len(samples), index=samples.index)
    for factor_name, factor_values in factors.items():
        if isinstance(factor_values, str):
            factor_values = [factor_values]
        elif not isinstance(factor_values, (list, tuple)):
            raise ValueError(
                f"Factor values for {factor_name} must be a string or a list/tuple."  # noqa: E501
            )
        mask &= samples[factor_name].isin(factor_values)
    matching_samples = samples[mask]["name"].astype(str).tolist()
    return matching_samples
