"""
Unit tests for allele harmonization logic.
"""

import pandas as pd
import pytest


def test_effect_allele_flip():
    """Harmonization must flip effect sign when alleles are swapped."""
    from bta.pgs.harmonize import flip_alleles

    row = {"effect_allele": "A", "other_allele": "T", "effect_weight": 0.5}
    flipped = flip_alleles(row, ref_effect="T", ref_other="A")
    assert flipped["effect_allele"] == "T"
    assert flipped["other_allele"] == "A"
    assert abs(flipped["effect_weight"] - (-0.5)) < 1e-9


def test_ambiguous_snps_removed():
    """A/T and C/G SNPs must be removed from main analysis."""
    from bta.pgs.harmonize import remove_ambiguous

    df = pd.DataFrame(
        {
            "effect_allele": ["A", "C", "G", "T", "A"],
            "other_allele": ["T", "G", "C", "A", "C"],
            "effect_weight": [0.1, 0.2, 0.3, 0.4, 0.5],
        }
    )
    clean = remove_ambiguous(df)
    assert len(clean) == 1
    assert clean.iloc[0]["effect_allele"] == "A"
    assert clean.iloc[0]["other_allele"] == "C"
