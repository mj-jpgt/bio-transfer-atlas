"""
Population score-shift metrics.
"""

import numpy as np
import pandas as pd


def standardized_shift(
    scores: pd.Series,
    group_labels: pd.Series,
    reference: str | None = None,
) -> pd.DataFrame:
    """
    Compute standardized population score shift.

    standardized_shift(p, ref) = |mean_p - mean_ref| / pooled_SD

    Parameters
    ----------
    scores : pd.Series
        Individual-level PGS scores.
    group_labels : pd.Series
        Population/group label per individual (same index as scores).
    reference : str, optional
        Reference population label. If None, uses global mean.

    Returns
    -------
    pd.DataFrame with columns [population, mean_score, sd_score, standardized_shift]
    """
    pooled_sd = scores.std(ddof=1)
    if pooled_sd == 0:
        pooled_sd = 1.0

    if reference is not None:
        ref_mean = scores[group_labels == reference].mean()
    else:
        ref_mean = scores.mean()

    rows = []
    for pop, grp in scores.groupby(group_labels):
        mean_p = grp.mean()
        sd_p = grp.std(ddof=1)
        shift = abs(mean_p - ref_mean) / pooled_sd
        rows.append(
            {
                "population": pop,
                "n": len(grp),
                "mean_score": mean_p,
                "sd_score": sd_p,
                "standardized_shift": shift,
            }
        )
    return pd.DataFrame(rows)


def bootstrap_shifts(
    scores: pd.Series,
    group_labels: pd.Series,
    reference: str | None = None,
    n_boot: int = 1000,
    seed: int = 719,
) -> pd.DataFrame:
    """
    Bootstrap confidence intervals for standardized population shift.
    """
    rng = np.random.default_rng(seed)
    boot_rows: dict[str, list[float]] = {}

    for _ in range(n_boot):
        idx = rng.choice(len(scores), size=len(scores), replace=True)
        s_b = scores.iloc[idx].reset_index(drop=True)
        g_b = group_labels.iloc[idx].reset_index(drop=True)
        df_b = standardized_shift(s_b, g_b, reference=reference)
        for _, row in df_b.iterrows():
            boot_rows.setdefault(row["population"], []).append(row["standardized_shift"])

    ci_rows = []
    for pop, vals in boot_rows.items():
        arr = np.array(vals)
        ci_rows.append(
            {
                "population": pop,
                "shift_mean": arr.mean(),
                "ci_lower": np.percentile(arr, 2.5),
                "ci_upper": np.percentile(arr, 97.5),
            }
        )
    return pd.DataFrame(ci_rows)
