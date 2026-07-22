from pathlib import Path

import pandas as pd


def test_master_table_schema():
    p = Path("data/modeling/master_variant_table.parquet")
    assert p.exists()
    d = pd.read_parquet(p)
    required = {
        "variant_id",
        "trait",
        "I2",
        "risk_class",
        "split_trait",
        "split_variant",
    }
    assert required.issubset(d.columns)


def test_leakage_free_variant_split():
    p = Path("data/modeling/master_variant_table.parquet")
    assert p.exists()
    d = pd.read_parquet(p, columns=["variant_id", "split_variant"]).drop_duplicates()
    counts = d.groupby("variant_id")["split_variant"].nunique()
    assert (counts == 1).all(), "variant leakage across split_variant folds"


def test_multisource_labels_schema_if_present():
    p = Path("data/labels/gwas_concordance_labels_multisource.parquet")
    if not p.exists():
        return
    d = pd.read_parquet(p)
    required = {"variant_id", "trait", "n_sources", "I2", "associated", "source_group"}
    assert required.issubset(d.columns)
