"""
Phase A2/A3: Evaluate AF_LD_SEL under LD-block holdout and reviewer baselines.

Baselines:
  - POP_DISTANCE: pairwise FST_*_EUR columns (Bitarello-style pop distance proxies)
  - RG_PROXY: block-level mean |AF_diff| and FST_like as cheap genetic-distance / discordance proxy
    (full Popcorn/LDSC can replace this when available)
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
        "--ld-blocks",
        default=str(ROOT / "data/modeling/ld_block_assignments_genomewide.parquet"),
    )
    p.add_argument(
        "--out",
        default=str(ROOT / "results/tables/ablation_ldblock_and_baselines_genomewide.csv"),
    )
    p.add_argument(
        "--max-train",
        type=int,
        default=250_000,
        help="Cap train rows to avoid HistGB MemoryError.",
    )
    p.add_argument("--max-sample", type=int, default=500_000)
    return p.parse_args()


def make_model() -> Pipeline:
    return Pipeline(
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


def cap_train(tr: pd.DataFrame, max_train: int) -> pd.DataFrame:
    if max_train <= 0 or len(tr) <= max_train:
        return tr
    y = tr["y_high_I2"].astype(int)
    pos = tr[y == 1]
    neg = tr[y == 0]
    n_pos = min(len(pos), max(1, int(max_train * (len(pos) / max(len(tr), 1)))))
    n_neg = min(len(neg), max_train - n_pos)
    parts = []
    if n_pos > 0:
        parts.append(pos.sample(n_pos, random_state=SEED))
    if n_neg > 0:
        parts.append(neg.sample(n_neg, random_state=SEED))
    return pd.concat(parts, ignore_index=True) if parts else tr.sample(max_train, random_state=SEED)


def eval_feats(
    df: pd.DataFrame, split_col: str, feats: list[str], name: str, max_train: int = 350_000
) -> dict:
    import gc

    use = [f for f in feats if f in df.columns]
    if not use:
        return {
            "split": split_col,
            "feature_group": name,
            "AUROC": np.nan,
            "AUPRC": np.nan,
            "n_feats": 0,
            "n_train": 0,
            "n_test": 0,
        }
    tr = cap_train(df[df[split_col] == "train"], max_train)
    te = df[df[split_col] == "test"]
    ytr = tr["y_high_I2"].astype(int).to_numpy()
    yte = te["y_high_I2"].astype(int).to_numpy()
    if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
        return {
            "split": split_col,
            "feature_group": name,
            "AUROC": np.nan,
            "AUPRC": np.nan,
            "n_feats": len(use),
            "n_train": len(tr),
            "n_test": len(te),
        }
    Xtr = tr[use].to_numpy(dtype=np.float32, copy=True)
    Xte = te[use].to_numpy(dtype=np.float32, copy=True)
    model = make_model()
    model.fit(Xtr, ytr)
    del Xtr
    gc.collect()
    p = model.predict_proba(Xte)[:, 1]
    del Xte, model
    gc.collect()
    return {
        "split": split_col,
        "feature_group": name,
        "AUROC": float(roc_auc_score(yte, p)),
        "AUPRC": float(average_precision_score(yte, p)),
        "n_feats": len(use),
        "n_train": len(tr),
        "n_test": len(te),
        "pos_rate_test": float(yte.mean()),
    }


def main() -> None:
    args = parse_args()
    groups = json.loads(Path(args.groups).read_text(encoding="utf-8"))
    af_ld_sel = groups["AF_LD_SEL"]
    fst = groups.get("FST", ["FST_like"])

    pop_distance = [
        c
        for c in [
            "FST_AFR_EUR",
            "FST_EAS_EUR",
            "FST_AMR_EUR",
            "FST_SAS_EUR",
            "FST_like",
            "AF_diff_AFR_EUR",
            "AF_diff_EAS_EUR",
            "AF_diff_AMR_EUR",
            "AF_diff_SAS_EUR",
        ]
        if True
    ]
    rg_proxy = [
        "FST_like",
        "AF_max_diff",
        "AF_var",
        "AIS",
        "LD_max_diff",
        "LD_entropy",
        "LD_EUR_AFR_ratio",
        "LD_EUR_EAS_ratio",
    ]

    need = set(
        ["variant_id", "trait", "y_high_I2", "split_variant", "associated"]
        + af_ld_sel
        + fst
        + pop_distance
        + rg_proxy
    )
    sample_path = ROOT / "data/modeling/_tmp_ldblock_associated_sample.parquet"
    print("Streaming associated sample to disk ...", flush=True)
    stream_sample_associated(
        Path(args.master),
        sorted(need),
        sample_path,
        max_rows=args.max_sample,
        mhc_keep_all=False,
        seed=SEED,
    )
    print("Loading sample ...", flush=True)
    df = pd.read_parquet(sample_path)
    print(f"Associated sample: {len(df):,}", flush=True)

    blocks = Path(args.ld_blocks)
    if not blocks.exists():
        raise SystemExit(f"Missing {blocks}; run assign_ld_block_splits.py first")
    ld = pd.read_parquet(blocks, columns=["variant_id", "split_ld_block", "ld_block"])
    df = df.merge(ld, on="variant_id", how="inner")
    print(f"After LD-block join: {len(df):,}", flush=True)

    # Block-level rg proxy: mean FST_like within block as extra feature
    block_rg = df.groupby("ld_block")["FST_like"].transform("mean")
    df["rg_block_proxy"] = block_rg
    rg_proxy = rg_proxy + ["rg_block_proxy"]

    # Trait-constant Z-concordance is NOT a peer of variant-varying AF/LD.
    # Emit separately as trait_constant_baseline (not in primary feature-group contest).
    z_feats = ["z_concordance_EUR_AFR", "z_concordance_EUR_EAS", "rg_EUR_AFR", "rg_EUR_EAS"]
    z_path = ROOT / "data/features/baselines/rg_real_by_trait.parquet"
    z_joined = []
    if z_path.exists() and "trait" in df.columns:
        ztab = pd.read_parquet(z_path)
        # Prefer z_concordance_* columns; fall back to legacy rg_* names
        avail = [c for c in z_feats if c in ztab.columns]
        if avail:
            df = df.merge(ztab[["trait"] + avail], on="trait", how="left")
            z_joined = avail
            print(f"Joined trait-constant concordance from {z_path}: {avail}", flush=True)

    # Peer contest must never include trait-constant concordance columns
    peer_cols = sorted(set(af_ld_sel + fst + pop_distance + rg_proxy))
    assert not any(
        "z_concordance" in c or (c.startswith("rg_") and c != "rg_block_proxy") for c in peer_cols
    ), "trait-constant concordance leaked into variant peer columns"

    rows = []
    max_train = args.max_train
    for split in ["split_variant", "split_ld_block"]:
        print(f"[{split}] AF_LD_SEL ...", flush=True)
        rows.append(eval_feats(df, split, af_ld_sel, "AF_LD_SEL", max_train=max_train))
        print(f"[{split}] FST ...", flush=True)
        rows.append(eval_feats(df, split, fst, "FST", max_train=max_train))
        print(f"[{split}] POP_DISTANCE ...", flush=True)
        rows.append(eval_feats(df, split, pop_distance, "POP_DISTANCE", max_train=max_train))
        print(f"[{split}] RG_PROXY ...", flush=True)
        rows.append(eval_feats(df, split, rg_proxy, "RG_PROXY", max_train=max_train))
        # Diagnostic only: trait-constant baseline (not a fair AF/LD comparator)
        if z_joined and all(c in df.columns for c in z_joined):
            print(f"[{split}] TRAIT_CONSTANT_Z (diagnostic; not peer contest) ...", flush=True)
            rows.append(
                eval_feats(df, split, z_joined, "TRAIT_CONSTANT_Z_DIAGNOSTIC", max_train=max_train)
            )

    out = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False, float_format="%.4f")
    print(out.to_string(index=False))
    print(f"Saved {out_path}")

    # Primary peer contest excludes TRAIT_CONSTANT_Z_DIAGNOSTIC
    peer = out[~out["feature_group"].astype(str).str.contains("TRAIT_CONSTANT", na=False)]
    peer_path = out_path.parent / "ablation_ldblock_peer_contest.csv"
    peer.to_csv(peer_path, index=False, float_format="%.4f")
    print(f"Saved peer contest (no trait-constant) -> {peer_path}")

    groups_path = Path(args.groups)
    g = json.loads(groups_path.read_text(encoding="utf-8"))
    g["POP_DISTANCE"] = [c for c in pop_distance if c in set(df.columns) and c != "rg_block_proxy"]
    g["RG_PROXY"] = [c for c in rg_proxy if c in df.columns or c == "rg_block_proxy"]
    if z_joined:
        g["TRAIT_CONSTANT_Z_DIAGNOSTIC"] = z_joined
        g.pop("RG_REAL", None)
    groups_path.write_text(json.dumps(g, indent=2), encoding="utf-8")
    print(f"Updated {groups_path} (RG_REAL removed from peer groups)")


if __name__ == "__main__":
    main()
