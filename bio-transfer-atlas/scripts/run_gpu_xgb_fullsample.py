#!/usr/bin/env python3
"""Larger-sample CUDA XGBoost on genomewide master (keeps A100 busy during intervention)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, roc_auc_score
from xgboost import XGBClassifier

ROOT = Path(__file__).resolve().parents[1]
SEED = 719
MAX_ROWS = 1_200_000


def main() -> None:
    groups = json.loads((ROOT / "data/modeling/feature_groups_genomewide_genomewide.json").read_text())
    feats = groups.get("AF_LD_SEL", [])
    mp = ROOT / "data/modeling/master_variant_table_genomewide_genomewide.parquet"
    ld = pd.read_parquet(
        ROOT / "data/modeling/ld_block_assignments_genomewide.parquet",
        columns=["variant_id", "split_ld_block"],
    )
    need = list(dict.fromkeys(["variant_id", "y_high_I2", "split_variant", *feats]))
    print(f"reading parquet cols={len(need)}", flush=True)
    df = pd.read_parquet(mp, columns=need)
    df = df[df["y_high_I2"].notna()]
    print(f"associated={len(df):,}", flush=True)
    if len(df) > MAX_ROWS:
        df = df.sample(MAX_ROWS, random_state=SEED)
    df = df.merge(ld.drop_duplicates("variant_id"), on="variant_id", how="inner")
    print(f"merged={len(df):,}", flush=True)
    use = [f for f in feats if f in df.columns]
    rows = []
    for split in ["split_ld_block", "split_variant"]:
        tr = df[df[split] == "train"]
        te = df[df[split] == "test"]
        imp = SimpleImputer(strategy="median")
        Xtr = imp.fit_transform(tr[use]).astype(np.float32)
        Xte = imp.transform(te[use]).astype(np.float32)
        ytr = tr["y_high_I2"].astype(int).to_numpy()
        yte = te["y_high_I2"].astype(int).to_numpy()
        clf = XGBClassifier(
            n_estimators=600,
            max_depth=10,
            learning_rate=0.04,
            subsample=0.85,
            colsample_bytree=0.85,
            min_child_weight=5,
            tree_method="hist",
            device="cuda",
            eval_metric="auc",
            n_jobs=8,
            random_state=SEED,
        )
        print(f"fit {split} n_train={len(tr):,} n_test={len(te):,}", flush=True)
        clf.fit(Xtr, ytr)
        p = clf.predict_proba(Xte)[:, 1]
        auroc = float(roc_auc_score(yte, p))
        auprc = float(average_precision_score(yte, p))
        print(f"[{split}] AUROC={auroc:.4f} AUPRC={auprc:.4f}", flush=True)
        rows.append(
            {
                "split": split,
                "feature_group": "AF_LD_SEL",
                "AUROC": auroc,
                "AUPRC": auprc,
                "n_train": len(tr),
                "n_test": len(te),
                "device": "cuda",
                "model": "xgb_hist_1p2M",
            }
        )
    out = ROOT / "results/tables/ablation_xgboost_gpu_fullsample.csv"
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.4f")
    print(f"Saved {out}", flush=True)
    print("GPU_FULLSAMPLE_DONE", flush=True)


if __name__ == "__main__":
    main()
