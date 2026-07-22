"""
Optional GPU XGBoost mirror of AF_LD_SEL / VEP_AF ablations (Lambda A100).
Primary paper model remains HistGB; this is a speed/GPU-utilization companion.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline

ROOT = Path(__file__).resolve().parents[1]
SEED = 719


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--sample",
        default=str(ROOT / "data/modeling/_tmp_ldblock_associated_sample.parquet"),
    )
    p.add_argument(
        "--ld-blocks",
        default=str(ROOT / "data/modeling/ld_block_assignments_genomewide.parquet"),
    )
    p.add_argument(
        "--groups",
        default=str(ROOT / "data/modeling/feature_groups_genomewide_genomewide.json"),
    )
    p.add_argument(
        "--out",
        default=str(ROOT / "results/tables/ablation_xgboost_gpu_companion.csv"),
    )
    p.add_argument("--device", default="cuda")
    return p.parse_args()


def main() -> None:
    try:
        from xgboost import XGBClassifier
    except ImportError as e:
        raise SystemExit(f"xgboost required: {e}")

    args = parse_args()
    groups = json.loads(Path(args.groups).read_text(encoding="utf-8"))
    df = pd.read_parquet(args.sample)
    ld = pd.read_parquet(args.ld_blocks, columns=["variant_id", "split_ld_block"])
    df = df.merge(ld, on="variant_id", how="inner")

    rg_path = ROOT / "data/features/baselines/rg_real_by_trait.parquet"
    if rg_path.exists() and "trait" in df.columns:
        rg = pd.read_parquet(rg_path)
        df = df.merge(rg, on="trait", how="left")

    feat_sets = {
        "AF_LD_SEL": groups.get("AF_LD_SEL", []),
        "VEP_AF": groups.get("VEP_AF", []),
    }
    # combine AF_LD_SEL+VEP if columns present
    if feat_sets["VEP_AF"]:
        feat_sets["AF_LD_SEL+VEP_AF"] = list(dict.fromkeys(feat_sets["AF_LD_SEL"] + feat_sets["VEP_AF"]))

    rows = []
    for split in ["split_ld_block"]:
        if split not in df.columns:
            continue
        tr = df[df[split] == "train"]
        te = df[df[split] == "test"]
        ytr = tr["y_high_I2"].astype(int).to_numpy()
        yte = te["y_high_I2"].astype(int).to_numpy()
        for name, feats in feat_sets.items():
            use = [f for f in feats if f in df.columns]
            if len(use) < 1 or len(np.unique(ytr)) < 2:
                continue
            tree_method = "hist"
            device = args.device
            clf = XGBClassifier(
                n_estimators=200,
                max_depth=6,
                learning_rate=0.08,
                subsample=0.8,
                colsample_bytree=0.8,
                random_state=SEED,
                tree_method=tree_method,
                device=device,
                eval_metric="auc",
            )
            pipe = Pipeline([("imputer", SimpleImputer(strategy="median")), ("xgb", clf)])
            try:
                pipe.fit(tr[use], ytr)
                p = pipe.predict_proba(te[use])[:, 1]
                auroc = float(roc_auc_score(yte, p))
            except Exception as exc:
                # CPU fallback
                clf.set_params(device="cpu")
                pipe = Pipeline([("imputer", SimpleImputer(strategy="median")), ("xgb", clf)])
                pipe.fit(tr[use], ytr)
                p = pipe.predict_proba(te[use])[:, 1]
                auroc = float(roc_auc_score(yte, p))
                device = f"cpu_fallback:{exc}"
            rows.append(
                {
                    "split": split,
                    "feature_group": name,
                    "AUROC": auroc,
                    "n_feats": len(use),
                    "n_train": len(tr),
                    "n_test": len(te),
                    "device": str(device),
                }
            )
            print(f"{split} {name}: AUROC={auroc:.4f} device={device}", flush=True)

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, float_format="%.4f")
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
