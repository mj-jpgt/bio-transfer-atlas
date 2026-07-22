#!/usr/bin/env python3
"""Held-out trait ranking: train AF_LD_SEL on some traits, evaluate heterogeneity ranking on held-out trait."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
SEED = 719


def main() -> None:
    sample = ROOT / "data/modeling/_tmp_ldblock_associated_sample.parquet"
    groups = json.loads(
        (ROOT / "data/modeling/feature_groups_genomewide_genomewide.json").read_text(encoding="utf-8")
    )
    out = ROOT / "results/tables/external_heldout_trait_ranking.csv"
    if not sample.exists():
        pd.DataFrame([{"status": "missing_sample"}]).to_csv(out, index=False)
        print("missing sample")
        return
    df = pd.read_parquet(sample)
    labels = ROOT / "data/labels/_tmp_associated_labels.parquet"
    if "trait" not in df.columns and labels.exists():
        lab = pd.read_parquet(labels, columns=["variant_id", "trait", "y_high_I2"])
        keep = ["variant_id", "trait"]
        if "y_high_I2" not in df.columns:
            keep.append("y_high_I2")
        df = df.merge(lab[keep].drop_duplicates("variant_id"), on="variant_id", how="inner")
    blocks = ROOT / "data/modeling/ld_block_assignments_genomewide.parquet"
    if "ld_block" not in df.columns and blocks.exists():
        ld = pd.read_parquet(blocks, columns=["variant_id", "ld_block"])
        df = df.merge(ld.drop_duplicates("variant_id"), on="variant_id", how="left")
    if "trait" not in df.columns or "y_high_I2" not in df.columns:
        pd.DataFrame([{"status": "need_trait_and_y_high_I2"}]).to_csv(out, index=False)
        print("need trait and y_high_I2")
        return
    af = [f for f in groups.get("AF_LD_SEL", []) if f in df.columns]
    fst = [f for f in groups.get("FST", ["FST_like"]) if f in df.columns]
    # Baselines
    feats_sets = {
        "AF_LD_SEL": af,
        "FST": fst,
        "MAF_proxy": [c for c in df.columns if "AF_" in c or "maf" in c.lower()][:5],
        "LD_proxy": [c for c in df.columns if "LD_" in c or "ldscore" in c.lower()][:5],
    }
    traits = sorted(df["trait"].astype(str).unique())
    rows = []
    for hold in traits:
        tr = df[df["trait"] != hold]
        te = df[df["trait"] == hold]
        if len(te) < 500 or len(tr) < 2000:
            continue
        ytr = tr["y_high_I2"].astype(int)
        yte = te["y_high_I2"].astype(int)
        if yte.nunique() < 2:
            continue
        # Random baseline
        rng = np.random.default_rng(SEED)
        p_rand = rng.random(len(te))
        rows.append(
            {
                "heldout_trait": hold,
                "model": "random",
                "AUROC": float(roc_auc_score(yte, p_rand)),
                "AUPRC": float(average_precision_score(yte, p_rand)),
                "n_test": len(te),
            }
        )
        for name, feats in feats_sets.items():
            feats = [f for f in feats if f in df.columns]
            if len(feats) < 1:
                continue
            imp = SimpleImputer(strategy="median")
            Xtr = imp.fit_transform(tr[feats])
            Xte = imp.transform(te[feats])
            clf = HistGradientBoostingClassifier(
                max_depth=6, learning_rate=0.05, max_iter=150, random_state=SEED
            )
            clf.fit(Xtr, ytr)
            p = clf.predict_proba(Xte)[:, 1]
            rows.append(
                {
                    "heldout_trait": hold,
                    "model": name,
                    "AUROC": float(roc_auc_score(yte, p)),
                    "AUPRC": float(average_precision_score(yte, p)),
                    "n_test": len(te),
                    "n_feats": len(feats),
                }
            )
            print(hold, name, rows[-1]["AUROC"], flush=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.4f")
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
