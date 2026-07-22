"""
Lean fine-map tier eval: uses cached associated sample only (no master re-stream).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, roc_auc_score
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
        "--tiers",
        default=str(ROOT / "data/labels/finemap_tiers_genomewide_zlead.parquet"),
    )
    p.add_argument(
        "--groups",
        default=str(ROOT / "data/modeling/feature_groups_genomewide_genomewide.json"),
    )
    p.add_argument(
        "--out",
        default=str(ROOT / "results/tables/ablation_finemap_tiers_zlead.csv"),
    )
    p.add_argument("--max-train", type=int, default=200_000)
    return p.parse_args()


def cap_train(tr: pd.DataFrame, max_train: int) -> pd.DataFrame:
    if len(tr) <= max_train:
        return tr
    y = tr["y_high_I2"].astype(int)
    pos, neg = tr[y == 1], tr[y == 0]
    n_pos = min(len(pos), max(1, int(max_train * len(pos) / max(len(tr), 1))))
    n_neg = min(len(neg), max_train - n_pos)
    parts = []
    if n_pos:
        parts.append(pos.sample(n_pos, random_state=SEED))
    if n_neg:
        parts.append(neg.sample(n_neg, random_state=SEED))
    return pd.concat(parts, ignore_index=True)


def eval_one(df: pd.DataFrame, feats: list[str], max_train: int) -> dict:
    tr = cap_train(df[df["split_variant"] == "train"], max_train)
    te = df[df["split_variant"] == "test"]
    ytr = tr["y_high_I2"].astype(int).to_numpy()
    yte = te["y_high_I2"].astype(int).to_numpy()
    if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2 or len(te) < 50:
        return {"AUROC": np.nan, "AUPRC": np.nan, "n_train": len(tr), "n_test": len(te)}
    use = [f for f in feats if f in df.columns]
    model = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            (
                "hgb",
                HistGradientBoostingClassifier(
                    random_state=SEED,
                    max_iter=80,
                    learning_rate=0.1,
                    max_depth=5,
                    max_bins=64,
                    early_stopping=True,
                    validation_fraction=0.1,
                    n_iter_no_change=8,
                ),
            ),
        ]
    )
    model.fit(tr[use].to_numpy(dtype=np.float32), ytr)
    p = model.predict_proba(te[use].to_numpy(dtype=np.float32))[:, 1]
    return {
        "AUROC": float(roc_auc_score(yte, p)),
        "AUPRC": float(average_precision_score(yte, p)),
        "n_train": len(tr),
        "n_test": len(te),
        "pos_rate_test": float(yte.mean()),
    }


def main() -> None:
    args = parse_args()
    sample = Path(args.sample)
    if not sample.exists():
        raise SystemExit(f"Missing cached sample {sample}")
    groups = json.loads(Path(args.groups).read_text(encoding="utf-8"))
    feats = groups["AF_LD_SEL"]
    print(f"Loading cached sample {sample} ...", flush=True)
    df = pd.read_parquet(sample)
    tiers = pd.read_parquet(args.tiers, columns=["variant_id", "trait", "finemap_tier"])
    if "trait" not in df.columns:
        # sample may be variant-level without trait — join on variant_id only (first trait)
        tiers = tiers.drop_duplicates("variant_id")
        df = df.merge(tiers, on="variant_id", how="inner")
    else:
        df = df.merge(tiers, on=["variant_id", "trait"], how="inner")
    print(f"Joined {len(df):,}", flush=True)

    rows = []
    for tier in ["fine_mapped", "tag_only", "all"]:
        sub = df if tier == "all" else df[df["finemap_tier"] == tier]
        print(f"tier={tier} n={len(sub):,} pos={sub['y_high_I2'].mean():.3f}", flush=True)
        m = eval_one(sub, feats, args.max_train)
        rows.append({"tier": tier, "feature_group": "AF_LD_SEL", "tier_method": "z_meta_lead", **m})

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, float_format="%.4f")
    print(out.to_string(index=False))
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
