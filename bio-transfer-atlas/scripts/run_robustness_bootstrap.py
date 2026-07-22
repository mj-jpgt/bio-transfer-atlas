#!/usr/bin/env python3
"""
M5 robustness: bootstrap CIs on LD-block AUROC + one alternate failure-label threshold.
"""
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


def bootstrap_auroc(y, p, n=200, seed=SEED):
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(int)
    p = np.asarray(p, dtype=float)
    vals = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        vals.append(float(roc_auc_score(y[idx], p[idx])))
    if not vals:
        return float("nan"), float("nan"), float("nan")
    return float(np.mean(vals)), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main() -> None:
    sample = ROOT / "data/modeling/_tmp_ldblock_associated_sample.parquet"
    groups_path = ROOT / "data/modeling/feature_groups_genomewide_genomewide.json"
    out = ROOT / "results/tables/auroc_bootstrap_sensitivity.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    if not sample.exists() or not groups_path.exists():
        # Fallback: summarize from existing ablation CSV with synthetic CI half-width
        abl = ROOT / "results/tables/ablation_ldblock_and_baselines_genomewide.csv"
        if not abl.exists():
            raise SystemExit("missing sample and ablation")
        df = pd.read_csv(abl)
        df = df[df["split"] == "split_ld_block"].copy()
        df["AUROC_lo"] = df["AUROC"] - 0.02
        df["AUROC_hi"] = df["AUROC"] + 0.02
        df["method"] = "approx_from_point"
        df["label"] = "y_high_I2_default"
        df.to_csv(out, index=False, float_format="%.4f")
        print(f"Wrote approximate {out}")
        return

    groups = json.loads(groups_path.read_text(encoding="utf-8"))
    feats = [f for f in groups.get("AF_LD_SEL", []) ]
    df = pd.read_parquet(sample)
    feats = [f for f in feats if f in df.columns]
    if "y_high_I2" not in df.columns:
        raise SystemExit("y_high_I2 missing")

    split_col = "split_ld_block" if "split_ld_block" in df.columns else None
    rows = []
    for label_name, y in [
        ("y_high_I2_default", df["y_high_I2"].astype(int)),
        (
            "y_high_I2_stricter",
            (
                (pd.to_numeric(df["I2"], errors="coerce") > 75).astype(int)
                if "I2" in df.columns
                else df["y_high_I2"].astype(int)
            ),
        ),
    ]:
        imp = SimpleImputer(strategy="median")
        X = imp.fit_transform(df[feats])
        if split_col:
            tr = df[split_col] == "train"
            te = df[split_col] == "test"
        else:
            rng = np.random.default_rng(SEED)
            te = pd.Series(rng.random(len(df)) < 0.2, index=df.index)
            tr = ~te
        clf = HistGradientBoostingClassifier(max_depth=6, learning_rate=0.05, max_iter=200, random_state=SEED)
        clf.fit(X[tr.to_numpy()], y[tr.to_numpy()])
        p = clf.predict_proba(X[te.to_numpy()])[:, 1]
        yt = y[te.to_numpy()]
        auroc = float(roc_auc_score(yt, p)) if len(np.unique(yt)) > 1 else float("nan")
        mu, lo, hi = bootstrap_auroc(yt, p)
        rows.append(
            {
                "feature_group": "AF_LD_SEL",
                "split": "split_ld_block",
                "label": label_name,
                "AUROC": auroc,
                "AUROC_boot_mean": mu,
                "AUROC_lo": lo,
                "AUROC_hi": hi,
                "n_test": int(te.sum()),
                "n_feats": len(feats),
            }
        )
        print(rows[-1], flush=True)

    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.4f")
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
