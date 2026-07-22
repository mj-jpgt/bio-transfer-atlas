"""
Phase B1: Evaluate AF_LD_SEL separately on fine_mapped vs tag_only tiers.
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

import sys

sys.path.insert(0, str(ROOT / "scripts"))
from _stream_sample_associated import stream_sample_associated  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--master",
        default=str(ROOT / "data/modeling/master_variant_table_genomewide_genomewide.parquet"),
    )
    p.add_argument(
        "--groups",
        default=str(ROOT / "data/modeling/feature_groups_genomewide_genomewide.json"),
    )
    p.add_argument(
        "--tiers",
        default=str(ROOT / "data/labels/finemap_tiers_genomewide.parquet"),
    )
    p.add_argument(
        "--out",
        default=str(ROOT / "results/tables/ablation_finemap_tiers.csv"),
    )
    p.add_argument("--max-train", type=int, default=250_000)
    p.add_argument("--max-sample", type=int, default=500_000)
    return p.parse_args()


def cap_train(tr: pd.DataFrame, max_train: int) -> pd.DataFrame:
    if max_train <= 0 or len(tr) <= max_train:
        return tr
    y = tr["y_high_I2"].astype(int)
    pos, neg = tr[y == 1], tr[y == 0]
    n_pos = min(len(pos), max(1, int(max_train * (len(pos) / max(len(tr), 1)))))
    n_neg = min(len(neg), max_train - n_pos)
    parts = []
    if n_pos:
        parts.append(pos.sample(n_pos, random_state=SEED))
    if n_neg:
        parts.append(neg.sample(n_neg, random_state=SEED))
    return pd.concat(parts, ignore_index=True) if parts else tr.sample(max_train, random_state=SEED)


def eval_one(df: pd.DataFrame, feats: list[str], max_train: int = 350_000) -> dict:
    import gc

    tr = cap_train(df[df["split_variant"] == "train"], max_train)
    te = df[df["split_variant"] == "test"]
    ytr = tr["y_high_I2"].astype(int).to_numpy()
    yte = te["y_high_I2"].astype(int).to_numpy()
    if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2 or len(te) < 50:
        return {
            "AUROC": np.nan,
            "AUPRC": np.nan,
            "n_train": len(tr),
            "n_test": len(te),
        }
    model = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            (
                "hgb",
                HistGradientBoostingClassifier(
                    random_state=SEED,
                    max_iter=120,
                    learning_rate=0.08,
                    max_depth=6,
                    max_bins=64,
                    early_stopping=True,
                    validation_fraction=0.1,
                    n_iter_no_change=10,
                ),
            ),
        ]
    )
    use = [f for f in feats if f in df.columns]
    Xtr = tr[use].to_numpy(dtype=np.float32, copy=True)
    Xte = te[use].to_numpy(dtype=np.float32, copy=True)
    model.fit(Xtr, ytr)
    del Xtr
    gc.collect()
    p = model.predict_proba(Xte)[:, 1]
    del Xte, model
    gc.collect()
    return {
        "AUROC": float(roc_auc_score(yte, p)),
        "AUPRC": float(average_precision_score(yte, p)),
        "n_train": len(tr),
        "n_test": len(te),
        "pos_rate_test": float(yte.mean()),
    }


def main() -> None:
    args = parse_args()
    groups = json.loads(Path(args.groups).read_text(encoding="utf-8"))
    feats = groups["AF_LD_SEL"]
    cols = ["variant_id", "trait", "y_high_I2", "split_variant", "associated"] + feats
    sample_path = ROOT / "data/modeling/_tmp_finemap_eval_sample.parquet"
    print("Streaming master sample ...", flush=True)
    stream_sample_associated(
        Path(args.master), cols, sample_path, max_rows=args.max_sample, seed=SEED
    )
    master = pd.read_parquet(sample_path)
    print(f"Sampled {len(master):,}", flush=True)
    tiers = pd.read_parquet(
        args.tiers, columns=["variant_id", "trait", "finemap_tier", "tier_method"]
    )
    df = master.merge(tiers, on=["variant_id", "trait"], how="inner")
    print(f"Joined {len(df):,}", flush=True)

    rows = []
    for tier in ["fine_mapped", "tag_only", "all"]:
        sub = df if tier == "all" else df[df["finemap_tier"] == tier]
        print(f"Eval tier={tier} n={len(sub):,}", flush=True)
        m = eval_one(sub, feats, max_train=args.max_train)
        rows.append({"tier": tier, "feature_group": "AF_LD_SEL", **m})

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, float_format="%.4f")
    print(out.to_string(index=False))
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
