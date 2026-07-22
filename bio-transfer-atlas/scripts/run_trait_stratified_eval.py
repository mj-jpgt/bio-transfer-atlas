"""
Science deepen: within-trait and leave-one-trait-out AF_LD_SEL evaluation.

Outputs:
  results/tables/ablation_per_trait_genomewide.csv
  results/tables/ablation_trait_holdout_genomewide.csv
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
RNG = np.random.default_rng(SEED)
TRAITS = ["CAD", "T2D", "BMI", "LDL"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Trait-stratified and leave-one-trait-out eval.")
    p.add_argument(
        "--master",
        default=str(ROOT / "data/modeling/master_variant_table_genomewide_genomewide.parquet"),
    )
    p.add_argument(
        "--groups",
        default=str(ROOT / "data/modeling/feature_groups_genomewide_genomewide.json"),
    )
    p.add_argument("--feature-group", default="AF_LD_SEL")
    p.add_argument("--tag", default="genomewide")
    p.add_argument("--n-boot", type=int, default=40)
    return p.parse_args()


def load_associated(master_path: Path, cols: list[str]) -> pd.DataFrame:
    import pyarrow.dataset as ds

    dataset = ds.dataset(str(master_path), format="parquet")
    available = set(dataset.schema.names)
    use = [c for c in cols if c in available]
    filt = ds.field("associated") == True  # noqa: E712
    scanner = dataset.scanner(columns=use, filter=filt, batch_size=500_000)
    chunks: list[pd.DataFrame] = []
    n = 0
    for batch in scanner.to_batches():
        chunks.append(batch.to_pandas())
        n += len(chunks[-1])
        if n % 2_000_000 < 500_000:
            print(f"  loaded {n:,} associated rows ...", flush=True)
    if not chunks:
        return pd.DataFrame()
    return pd.concat(chunks, ignore_index=True)


def make_model() -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            (
                "hgb",
                HistGradientBoostingClassifier(
                    random_state=SEED,
                    max_iter=250,
                    learning_rate=0.05,
                ),
            ),
        ]
    )


def bootstrap_auroc(y: np.ndarray, p: np.ndarray, n_boot: int) -> tuple[float, float]:
    vals = []
    n = len(y)
    for _ in range(n_boot):
        idx = RNG.integers(0, n, n)
        yy = y[idx]
        if len(np.unique(yy)) < 2:
            continue
        vals.append(roc_auc_score(yy, p[idx]))
    if not vals:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def eval_split(
    tr: pd.DataFrame,
    te: pd.DataFrame,
    feats: list[str],
    n_boot: int,
    permute_train: bool = False,
) -> dict:
    ytr = tr["y_high_I2"].astype(int).to_numpy()
    yte = te["y_high_I2"].astype(int).to_numpy()
    if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
        return {
            "AUROC": np.nan,
            "AUROC_lo": np.nan,
            "AUROC_hi": np.nan,
            "AUPRC": np.nan,
            "n_train": len(tr),
            "n_test": len(te),
            "pos_rate_test": float(yte.mean()) if len(yte) else np.nan,
        }
    Xtr = tr[feats].to_numpy(dtype=np.float32)
    Xte = te[feats].to_numpy(dtype=np.float32)
    yfit = RNG.permutation(ytr) if permute_train else ytr
    model = make_model()
    model.fit(Xtr, yfit)
    p = model.predict_proba(Xte)[:, 1]
    au = float(roc_auc_score(yte, p))
    ap = float(average_precision_score(yte, p))
    lo, hi = bootstrap_auroc(yte, p, n_boot)
    return {
        "AUROC": au,
        "AUROC_lo": lo,
        "AUROC_hi": hi,
        "AUPRC": ap,
        "n_train": len(tr),
        "n_test": len(te),
        "pos_rate_test": float(yte.mean()),
    }


def main() -> None:
    args = parse_args()
    groups = json.loads(Path(args.groups).read_text(encoding="utf-8"))
    feats = groups[args.feature_group]
    cols = ["trait", "y_high_I2", "split_variant", "associated"] + feats
    print(f"Loading associated rows from {args.master} ...", flush=True)
    df = load_associated(Path(args.master), cols)
    print(f"Associated rows: {len(df):,}", flush=True)
    if "trait" not in df.columns:
        raise SystemExit("master missing trait column")

    out_dir = ROOT / "results/tables"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Within-trait, split_variant ---
    per_trait_rows = []
    for trait in TRAITS:
        sub = df[df["trait"] == trait]
        tr = sub[sub["split_variant"] == "train"]
        te = sub[sub["split_variant"] == "test"]
        print(f"[within-trait] {trait} train={len(tr):,} test={len(te):,}", flush=True)
        metrics = eval_split(tr, te, feats, args.n_boot, permute_train=False)
        per_trait_rows.append(
            {
                "experiment": "within_trait",
                "trait": trait,
                "feature_group": args.feature_group,
                "split": "split_variant",
                "permuted": False,
                **metrics,
            }
        )
        metrics_p = eval_split(tr, te, feats, args.n_boot, permute_train=True)
        per_trait_rows.append(
            {
                "experiment": "within_trait",
                "trait": trait,
                "feature_group": "PERMUTED",
                "split": "split_variant",
                "permuted": True,
                **metrics_p,
            }
        )

    per_trait = pd.DataFrame(per_trait_rows)
    per_path = out_dir / f"ablation_per_trait_{args.tag}.csv"
    per_trait.to_csv(per_path, index=False, float_format="%.4f")
    real = per_trait[(~per_trait["permuted"])]
    print(
        f"Within-trait AF_LD_SEL AUROC mean={real['AUROC'].mean():.3f} "
        f"min={real['AUROC'].min():.3f} max={real['AUROC'].max():.3f}"
    )
    print(f"Saved {per_path}")

    # --- Leave-one-trait-out ---
    holdout_rows = []
    for held in TRAITS:
        tr = df[df["trait"] != held]
        te = df[df["trait"] == held]
        print(f"[LOTO] holdout={held} train={len(tr):,} test={len(te):,}", flush=True)
        metrics = eval_split(tr, te, feats, args.n_boot, permute_train=False)
        holdout_rows.append(
            {
                "experiment": "leave_one_trait_out",
                "held_out_trait": held,
                "feature_group": args.feature_group,
                "permuted": False,
                **metrics,
            }
        )
        metrics_p = eval_split(tr, te, feats, args.n_boot, permute_train=True)
        holdout_rows.append(
            {
                "experiment": "leave_one_trait_out",
                "held_out_trait": held,
                "feature_group": "PERMUTED",
                "permuted": True,
                **metrics_p,
            }
        )

    holdout = pd.DataFrame(holdout_rows)
    ho_path = out_dir / f"ablation_trait_holdout_{args.tag}.csv"
    holdout.to_csv(ho_path, index=False, float_format="%.4f")
    real_ho = holdout[~holdout["permuted"]]
    print(
        f"Trait-holdout AF_LD_SEL AUROC mean={real_ho['AUROC'].mean():.3f} "
        f"min={real_ho['AUROC'].min():.3f} max={real_ho['AUROC'].max():.3f}"
    )
    print(f"Saved {ho_path}")


if __name__ == "__main__":
    main()
