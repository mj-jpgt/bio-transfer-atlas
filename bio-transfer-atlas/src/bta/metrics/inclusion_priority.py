"""
Inclusion Priority Score — composite index over ancestry space.

InclusionPriority =
    rank_norm(genetic_distance)
  + rank_norm(rank_instability)
  + rank_norm(pathway_instability)
  + rank_norm(variant_missingness)
  + rank_norm(low_density)
"""

from __future__ import annotations

import pandas as pd


def rank_norm(series: pd.Series) -> pd.Series:
    return series.rank() / len(series)


def compute_inclusion_priority(
    df: pd.DataFrame,
    components: list[str] | None = None,
    weights: dict[str, float] | None = None,
) -> pd.Series:
    """
    Compute composite inclusion priority score.

    Parameters
    ----------
    df         : DataFrame with component columns.
    components : Column names to include (default: all numeric).
    weights    : Optional per-component weights (default: equal).

    Returns
    -------
    pd.Series of priority scores (higher = higher inclusion priority).
    """
    if components is None:
        components = df.select_dtypes("number").columns.tolist()

    if weights is None:
        weights = {c: 1.0 for c in components}

    priority = sum(weights[c] * rank_norm(df[c]) for c in components)
    return priority
