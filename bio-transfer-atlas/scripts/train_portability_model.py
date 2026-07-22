"""
Phase 18.1: Train and persist portability-risk model (HGB, AF_LD_SEL).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline

from intervention_common import GW_MASTER, GW_MODEL, GW_TAG, ROOT, SEED

MDIR = ROOT / "data/modeling"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Train portability model and save joblib.")
    p.add_argument("--master", default=str(GW_MASTER))
    p.add_argument("--groups", default=str(MDIR / f"feature_groups_genomewide_{GW_TAG}.json"))
    p.add_argument("--feature-group", default="AF_LD_SEL")
    p.add_argument("--target", default="y_high_I2")
    p.add_argument("--split", default="split_variant")
    p.add_argument("--out", default=str(GW_MODEL))
    return p.parse_args()


def load_associated(master_path: Path, cols: list[str]) -> pd.DataFrame:
    import pyarrow.dataset as ds

    dataset = ds.dataset(str(master_path), format="parquet")
    filt = ds.field("associated") == True  # noqa: E712
    scanner = dataset.scanner(columns=cols, filter=filt, batch_size=500_000)
    chunks = []
    n = 0
    for batch in scanner.to_batches():
        chunks.append(batch.to_pandas())
        n += len(chunks[-1])
        if n % 2_000_000 < 500_000:
            print(f"  loaded {n:,} associated rows ...", flush=True)
    if not chunks:
        return pd.DataFrame(columns=cols)
    return pd.concat(chunks, ignore_index=True)


def main() -> None:
    args = parse_args()
    out_path = Path(args.out)
    meta_path = out_path.with_suffix(".meta.json")
    groups_path = Path(args.groups)
    if not groups_path.exists():
        groups_path = MDIR / "feature_groups.chr22.json"
    groups = json.loads(groups_path.read_text())
    feats = groups[args.feature_group]
    print(f"Features ({args.feature_group}): {len(feats)}")

    master_path = Path(args.master)
    cols = [args.split, args.target] + feats
    master = load_associated(master_path, cols)
    print(f"Associated rows: {len(master):,}")

    train = master[master[args.split] == "train"]
    test = master[master[args.split] == "test"]
    print(f"Train: {len(train):,}  Test: {len(test):,}")

    Xtr = train[feats].to_numpy(dtype=np.float32)
    ytr = train[args.target].to_numpy()
    Xte = test[feats].to_numpy(dtype=np.float32)
    yte = test[args.target].to_numpy()

    model = make_pipeline(
        SimpleImputer(strategy="median"),
        HistGradientBoostingClassifier(
            random_state=SEED,
            max_iter=300,
            learning_rate=0.05,
        ),
    )
    model.fit(Xtr, ytr)
    proba = model.predict_proba(Xte)[:, 1]
    auroc = float(roc_auc_score(yte, proba))
    print(f"Test AUROC ({args.split}): {auroc:.4f}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, out_path)
    meta = {
        "feature_group": args.feature_group,
        "features": feats,
        "target": args.target,
        "split": args.split,
        "seed": SEED,
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "test_auroc": auroc,
        "master": str(args.master),
    }
    meta_path.write_text(json.dumps(meta, indent=2))
    print(f"Saved model -> {out_path}")
    print(f"Saved meta  -> {meta_path}")


if __name__ == "__main__":
    main()
