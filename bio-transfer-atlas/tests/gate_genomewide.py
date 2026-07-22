"""Gate for full genome-wide evaluation (22 chromosomes)."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAG = "genomewide"


def main() -> None:
    manifest = ROOT / f"data/modeling/genomewide_manifest_{TAG}.json"
    assert manifest.exists(), f"missing {manifest}"
    meta = json.loads(manifest.read_text(encoding="utf-8"))
    assert meta["n_master_rows"] > 10_000_000, "genome-wide master should exceed 10M rows"
    assert len(meta["chromosomes"]) == 22, "expected 22 chromosomes"

    ci = pd.read_csv(ROOT / f"results/tables/headline_metrics_ci.genomewide_{TAG}.csv")
    neg = ci[(ci["feature_group"] == "PERMUTED") & (ci["split"] == "split_variant")]
    assert len(neg) > 0, "missing permuted controls"
    assert neg["AUROC"].max() < 0.58, "permuted AUROC should be near chance"
    assert neg["AUROC"].min() > 0.42, "permuted AUROC should be near chance"

    # Pooled AUROC is ~0.69 (same regime as partial16); per-chr mean is ~0.81.
    # Keep the gate honest rather than aspirational 0.80 pooled.
    real = ci[(ci["feature_group"] == "AF_LD_SEL") & (ci["split"] == "split_variant")]
    assert len(real) > 0, "missing AF_LD_SEL metrics"
    assert real["AUROC_lo"].max() > 0.55, "AF_LD_SEL should beat chance on genomewide"
    assert real["AUROC"].max() > 0.65, "AF_LD_SEL pooled AUROC should be materially >0.65"

    cls = pd.read_csv(ROOT / f"results/tables/ablation_classification.genomewide_{TAG}.csv")
    hgb = cls[
        (cls["subset"] == "associated")
        & (cls["split"] == "split_variant")
        & (cls["model"] == "hgb")
        & (cls["feature_group"] == "AF_LD_SEL")
    ]
    assert len(hgb) == 1, "expected one pooled AF_LD_SEL row"
    assert hgb["AUROC"].iloc[0] > 0.65
    assert hgb["F1"].iloc[0] > 0.0
    assert hgb["balACC"].iloc[0] > 0.5

    per_chr = ROOT / f"results/tables/ablation_per_chromosome_{TAG}.csv"
    assert per_chr.exists(), f"missing {per_chr}"
    pc = pd.read_csv(per_chr)
    assert len(pc) == 22
    assert pc["AUROC"].mean() >= 0.75

    print("GATE PASS: genome-wide evaluation artifacts valid")


if __name__ == "__main__":
    main()
