"""Gate for science-deepen holdout / trait-stratified evaluation artifacts."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAG = "genomewide"


def main() -> None:
    per_trait = ROOT / f"results/tables/ablation_per_trait_{TAG}.csv"
    holdout = ROOT / f"results/tables/ablation_trait_holdout_{TAG}.csv"
    headline = ROOT / "results/tables/headline_metrics_ci.genomewide_holdouts.csv"
    pooled = ROOT / f"results/tables/headline_metrics_ci.genomewide_{TAG}.csv"

    assert per_trait.exists(), f"missing {per_trait}"
    assert holdout.exists(), f"missing {holdout}"
    assert headline.exists(), f"missing {headline}"

    pt = pd.read_csv(per_trait)
    if pt["permuted"].dtype == object:
        pt["permuted"] = pt["permuted"].astype(str).str.lower().isin(["true", "1"])
    real = pt[(pt["feature_group"] == "AF_LD_SEL") & (~pt["permuted"])]
    assert len(real) >= 4, "expected within-trait rows for 4 traits"
    assert real["AUROC"].notna().all()

    pooled_ci = pd.read_csv(pooled)
    pooled_au = pooled_ci[
        (pooled_ci["feature_group"] == "AF_LD_SEL")
        & (pooled_ci["split"] == "split_variant")
        & (pooled_ci["subset"] == "associated")
    ]["AUROC"].max()
    assert real["AUROC"].mean() > float(pooled_au), (
        f"within-trait mean AUROC {real['AUROC'].mean():.3f} should exceed pooled {pooled_au:.3f}"
    )

    ho = pd.read_csv(holdout)
    if ho["permuted"].dtype == object:
        ho["permuted"] = ho["permuted"].astype(str).str.lower().isin(["true", "1"])
    real_ho = ho[(ho["feature_group"] == "AF_LD_SEL") & (~ho["permuted"])]
    perm_ho = ho[ho["permuted"]]
    assert len(real_ho) >= 4, "expected LOTO rows for 4 traits"
    assert len(perm_ho) >= 4, "expected LOTO permutation controls"
    # Publishable either way: just require permutation control near chance on average
    assert perm_ho["AUROC"].mean() < 0.60, "LOTO permutation should be near chance"

    hl = pd.read_csv(headline)
    splits = set(hl["split"].unique())
    assert "split_trait" in splits, "holdout headline missing split_trait"
    assert "split_source" in splits or "split_variant" in splits

    print("GATE PASS: science-deepen holdout artifacts valid")


if __name__ == "__main__":
    main()
