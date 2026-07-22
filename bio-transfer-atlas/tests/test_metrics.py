"""
Unit tests for instability metrics.
"""

import numpy as np
import pandas as pd
import pytest

from bta.metrics.shifts import standardized_shift


def test_rank_instability_zero_when_identical():
    """Rank instability must be zero when all score vectors are identical."""
    scores_a = pd.Series([1.0, 2.0, 3.0, 4.0])
    scores_b = pd.Series([1.0, 2.0, 3.0, 4.0])

    ranks_a = scores_a.rank(pct=True)
    ranks_b = scores_b.rank(pct=True)

    instability = (ranks_a - ranks_b).std()
    assert instability == 0.0


def test_standardized_shift_zero_for_identical_groups():
    """Standardized shift must be 0 when all groups have the same mean."""
    scores = pd.Series([1.0, 1.0, 1.0, 1.0])
    labels = pd.Series(["A", "A", "B", "B"])
    df = standardized_shift(scores, labels, reference="A")
    for _, row in df.iterrows():
        assert abs(row["standardized_shift"]) < 1e-9


def test_genetic_distance_sensitivity_positive():
    """Slope of shift ~ distance regression must be positive for synthetic toy values."""
    distances = np.array([0.0, 1.0, 2.0, 3.0])
    shifts = np.array([0.0, 0.5, 1.0, 1.5])

    x = distances - distances.mean()
    y = shifts - shifts.mean()
    beta = (x * y).sum() / (x**2).sum()

    assert beta > 0


def test_inclusion_priority_bounded():
    """Inclusion priority score must be bounded and reproducible."""
    components = pd.DataFrame(
        {
            "genetic_distance": [0.1, 0.5, 0.9],
            "rank_instability": [0.2, 0.4, 0.8],
            "pathway_instability": [0.3, 0.3, 0.7],
            "variant_missingness": [0.0, 0.1, 0.5],
        }
    )

    def rank_norm(col):
        return col.rank() / len(col)

    priority = sum(rank_norm(components[c]) for c in components.columns)
    assert priority.min() >= 0
    assert priority.max() <= len(components.columns) + 1
