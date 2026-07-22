"""
Phase A1: MHC stratified sensitivity for AF_LD_SEL mechanism model.

MHC (GRCh38): chr6:28510020-33480577
Reports AUROC with MHC excluded, MHC-only, and full + is_mhc indicator.

Low-RAM path: stream-sample associated rows to a temp parquet, then fit capped HGB.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from _stream_sample_associated import stream_sample_associated  # noqa: E402

SEED = 719
MHC_CHR = "6"
MHC_START = 28_510_020
MHC_END = 33_480_577


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="MHC sensitivity for portability model.")
    p.add_argument(
        "--master",
        default=str(ROOT / "data/modeling/master_variant_table_genomewide_genomewide.parquet"),
    )
    p.add_argument(
        "--groups",
        default=str(ROOT / "data/modeling/feature_groups_genomewide_genomewide.json"),
    )
    p.add_argument("--feature-group", default="AF_LD_SEL")
    p.add_argument("--max-sample", type=int, default=450_000, help="Non-MHC associated rows to keep")
    p.add_argument("--max-train", type=int, default=250_000)
    p.add_argument("--out", default=str(ROOT / "results/tables/mhc_sensitivity_genomewide.csv"))
    return p.parse_args()


def parse_pos(vid: str) -> tuple[str, int] | None:
    parts = str(vid).split(":")
    if len(parts) < 2:
        return None
    try:
        return parts[0], int(parts[1])
    except ValueError:
        return None


def is_mhc(vid: str) -> bool:
    parsed = parse_pos(vid)
    if parsed is None:
        return False
    chrom, pos = parsed
    return chrom == MHC_CHR and MHC_START <= pos <= MHC_END


def make_model() -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            (
                "hgb",
                HistGradientBoostingClassifier(
                    random_state=SEED,
                    max_iter=100,
                    learning_rate=0.08,
                    max_depth=5,
                    max_bins=64,
                    early_stopping=True,
                    validation_fraction=0.1,
                    n_iter_no_change=8,
                ),
            ),
        ]
    )


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


def eval_split(
    tr: pd.DataFrame, te: pd.DataFrame, feats: list[str], max_train: int = 250_000
) -> dict:
    import gc

    tr = cap_train(tr, max_train)
    ytr = tr["y_high_I2"].astype(int).to_numpy()
    yte = te["y_high_I2"].astype(int).to_numpy()
    if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
        return {"AUROC": np.nan, "AUPRC": np.nan, "n_train": len(tr), "n_test": len(te)}
    Xtr = np.ascontiguousarray(tr[feats].to_numpy(dtype=np.float32))
    Xte = np.ascontiguousarray(te[feats].to_numpy(dtype=np.float32))
    model = make_model()
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
        "n_mhc_test": int(te["is_mhc"].sum()) if "is_mhc" in te.columns else 0,
    }


def main() -> None:
    args = parse_args()
    groups = json.loads(Path(args.groups).read_text(encoding="utf-8"))
    feats = list(groups[args.feature_group])
    cols = ["variant_id", "y_high_I2", "split_variant", "associated"] + feats

    sample_path = ROOT / "data/modeling/_tmp_mhc_associated_sample.parquet"
    if sample_path.exists() and sample_path.stat().st_size > 1_000_000:
        print(f"Reusing existing sample {sample_path}", flush=True)
    else:
        print("Streaming associated sample to disk ...", flush=True)
        stream_sample_associated(
            Path(args.master),
            cols,
            sample_path,
            max_rows=args.max_sample,
            mhc_keep_all=True,
            seed=SEED,
        )
    print("Loading sample parquet ...", flush=True)
    df = pd.read_parquet(sample_path)
    print(f"Associated sample: {len(df):,}", flush=True)
    df["is_mhc"] = df["variant_id"].map(is_mhc)
    n_mhc = int(df["is_mhc"].sum())
    print(f"MHC rows: {n_mhc:,} ({100 * n_mhc / max(len(df), 1):.2f}%)", flush=True)

    rows = []
    max_train = args.max_train

    print("Eval full_pooled ...", flush=True)
    m = eval_split(
        df[df["split_variant"] == "train"],
        df[df["split_variant"] == "test"],
        feats,
        max_train=max_train,
    )
    rows.append({"setting": "full_pooled", "feature_group": args.feature_group, **m})

    print("Eval mhc_excluded ...", flush=True)
    df_ex = df[~df["is_mhc"]]
    m = eval_split(
        df_ex[df_ex["split_variant"] == "train"],
        df_ex[df_ex["split_variant"] == "test"],
        feats,
        max_train=max_train,
    )
    rows.append({"setting": "mhc_excluded", "feature_group": args.feature_group, **m})
    del df_ex

    print("Eval mhc_only ...", flush=True)
    df_m = df[df["is_mhc"]]
    m = eval_split(
        df_m[df_m["split_variant"] == "train"],
        df_m[df_m["split_variant"] == "test"],
        feats,
        max_train=0,
    )
    rows.append({"setting": "mhc_only", "feature_group": args.feature_group, **m})
    del df_m

    print("Eval full_plus_is_mhc ...", flush=True)
    feats_ind = feats + ["is_mhc"]
    df["is_mhc"] = df["is_mhc"].astype(np.float32)
    m = eval_split(
        df[df["split_variant"] == "train"],
        df[df["split_variant"] == "test"],
        feats_ind,
        max_train=max_train,
    )
    rows.append(
        {"setting": "full_plus_is_mhc", "feature_group": args.feature_group + "+is_mhc", **m}
    )

    out = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False, float_format="%.4f")
    print(out.to_string(index=False))
    print(f"Saved {out_path}")
    # note sampling in sidecar
    note = out_path.with_suffix(".note.txt")
    note.write_text(
        f"Sampled non-MHC associated rows (max={args.max_sample}) + all MHC; "
        f"train capped at {args.max_train}. Not full 3M-row fit.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
