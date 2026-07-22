"""
Rank instability metric across multiple PGS versions for the same trait.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def compute_rank_instability(
    score_df: pd.DataFrame,
    score_cols: list[str],
    group_col: str | None = None,
) -> pd.DataFrame:
    """
    rank_instability_i = SD(percentile_rank_i across score versions)

    Parameters
    ----------
    score_df  : DataFrame with individuals as rows, score versions as columns.
    score_cols: List of column names for different score versions.
    group_col : Optional column for population label.

    Returns
    -------
    pd.DataFrame with 'rank_instability' per individual,
    plus optional population-level summary.
    """
    ranks = score_df[score_cols].rank(pct=True)
    instability = ranks.std(axis=1, ddof=1)

    out = score_df[[group_col]].copy() if group_col else pd.DataFrame(index=score_df.index)
    out["rank_instability"] = instability

    if group_col:
        summary = (
            out.groupby(group_col)["rank_instability"]
            .agg(
                mean_instability="mean",
                p95_instability=lambda x: np.percentile(x, 95),
                n="count",
            )
            .reset_index()
        )
        return out, summary

    return out
