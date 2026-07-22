#!/usr/bin/env python3
"""Grouped family permutation importance on LD-block holdout (primary attribution)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
SEED = 719
N_REPEATS = 50


def main() -> None:
    sample = ROOT / "data/modeling/_tmp_ldblock_associated_sample.parquet"
    groups = json.loads(
        (ROOT / "data/modeling/feature_groups_genomewide_genomewide.json").read_text(encoding="utf-8")
    )
    out = ROOT / "results/tables/grouped_permutation_importance_ldblock.csv"
    if not sample.exists():
        pd.DataFrame([{"status": "missing_sample"}]).to_csv(out, index=False)
        return
    df = pd.read_parquet(sample)
    blocks = ROOT / "data/modeling/ld_block_assignments_genomewide.parquet"
    if "split_ld_block" not in df.columns and blocks.exists():
        ld = pd.read_parquet(blocks, columns=["variant_id", "split_ld_block"])
        df = df.merge(ld.drop_duplicates("variant_id"), on="variant_id", how="inner")

    families = {
        "AF": [c for c in groups.get("AF", []) if c in df.columns],
        "LD": [c for c in groups.get("LD", []) if c in df.columns],
        "SEL": [c for c in groups.get("SEL", []) if c in df.columns],
        "distance": [c for c in groups.get("POP_DISTANCE", []) if c in df.columns],
        "FST": [c for c in groups.get("FST", []) if c in df.columns],
    }
    all_feats = sorted(set(sum(families.values(), [])))
    assert not any("z_concordance" in c or c.startswith("rg_") for c in all_feats)

    tr = df[df["split_ld_block"] == "train"]
    te = df[df["split_ld_block"] == "test"]
    ytr = tr["y_high_I2"].astype(int)
    yte = te["y_high_I2"].astype(int)
    imp = SimpleImputer(strategy="median")
    Xtr = imp.fit_transform(tr[all_feats])
    Xte = imp.transform(te[all_feats])
    clf = HistGradientBoostingClassifier(
        max_depth=6, learning_rate=0.05, max_iter=200, random_state=SEED
    )
    clf.fit(Xtr, ytr)
    base = float(roc_auc_score(yte, clf.predict_proba(Xte)[:, 1]))
    feat_index = {c: i for i, c in enumerate(all_feats)}
    rng = np.random.default_rng(SEED)
    rows = []
    for fam, cols in families.items():
        if not cols:
            continue
        idxs = [feat_index[c] for c in cols if c in feat_index]
        drops = []
        for _ in range(N_REPEATS):
            Xp = Xte.copy()
            for j in idxs:
                Xp[:, j] = rng.permutation(Xp[:, j])
            drops.append(base - float(roc_auc_score(yte, clf.predict_proba(Xp)[:, 1])))
        rows.append(
            {
                "feature_family": fam,
                "n_features": len(idxs),
                "baseline_auroc": base,
                "mean_auroc_drop": float(np.mean(drops)),
                "std_auroc_drop": float(np.std(drops)),
                "n_repeats": N_REPEATS,
                "attribution": "grouped_permutation",
                "note": "descriptive_attention_weights are not primary attribution",
            }
        )
        print(fam, rows[-1]["mean_auroc_drop"], flush=True)
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.6g")
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
