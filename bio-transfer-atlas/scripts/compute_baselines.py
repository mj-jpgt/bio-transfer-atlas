"""
FAIRGEN-Open Stage 11: Mechanism-ablation baseline models
"""
from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, average_precision_score, f1_score, balanced_accuracy_score,
    mean_squared_error, mean_absolute_error,
)
from scipy import stats

root    = Path(__file__).resolve().parents[1]
mdir    = root / "data/modeling"
out_dir = root / "results/tables"
out_dir.mkdir(parents=True, exist_ok=True)
SEED = 719


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--master", default=str(mdir / "master_variant_table.parquet"))
    p.add_argument("--groups", default=str(mdir / "feature_groups.json"))
    p.add_argument("--subsets", default="all,associated", help="Comma list: all, associated")
    p.add_argument("--associated-only", action="store_true", help="Load only associated rows (saves RAM)")
    p.add_argument("--splits", default="split_variant,split_trait", help="Comma list of split columns")
    p.add_argument("--skip-riskclass", action="store_true", help="Skip 3-class risk model (saves RAM)")
    p.add_argument("--out-suffix", default="", help="Suffix for output CSVs, e.g. .genomewide_partial6")
    return p.parse_args()


def needed_columns(groups: dict, splits: list[str], skip_riskclass: bool) -> list[str]:
    cols: list[str] = ["associated", "y_high_I2", "I2"]
    if not skip_riskclass:
        cols.append("y_risk_class")
    cols.extend(splits)
    for feats in groups.values():
        cols.extend(feats)
    # preserve order, drop dups
    seen: set[str] = set()
    out: list[str] = []
    for c in cols:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out


def load_master(
    path: str,
    associated_only: bool,
    columns: list[str] | None = None,
) -> pd.DataFrame:
    """Load master parquet; stream associated rows to avoid OOM on large tables."""
    if not associated_only:
        return pd.read_parquet(path, columns=columns)
    import pyarrow.dataset as ds

    dataset = ds.dataset(path, format="parquet")
    available = set(dataset.schema.names)
    cols = [c for c in (columns or list(available)) if c in available]
    filt = ds.field("associated") == True  # noqa: E712
    scanner = dataset.scanner(columns=cols, filter=filt, batch_size=500_000)
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


args = parse_args()
groups = json.loads(Path(args.groups).read_text())
split_cols_early = [s.strip() for s in args.splits.split(",") if s.strip()]
cols = needed_columns(groups, split_cols_early, args.skip_riskclass)
master = load_master(args.master, args.associated_only, columns=cols)
print(f"master: {len(master):,} rows from {args.master} ({len(cols)} cols)")


def ece(y, p, bins=10):
    edges = np.linspace(0, 1, bins + 1)
    e = 0.0
    for i in range(bins):
        m = (p >= edges[i]) & (p < edges[i + 1])
        if m.sum() == 0:
            continue
        e += (m.mean()) * abs(y[m].mean() - p[m].mean())
    return float(e)


def get_xy(df, feats, target):
    X = df[feats].to_numpy(dtype=np.float32)
    y = df[target].to_numpy()
    return X, y


def run_classification(target, split_col, subset_name, df):
    rows = []
    tr = df[df[split_col] == "train"]
    te = df[df[split_col] == "test"]
    for gname, feats in groups.items():
        if not feats:
            continue
        Xtr, ytr = get_xy(tr, feats, target)
        Xte, yte = get_xy(te, feats, target)
        if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
            continue
        for mname, model in [
            ("logreg", make_pipeline(
                SimpleImputer(strategy="median"), StandardScaler(),
                LogisticRegression(max_iter=2000, class_weight="balanced", C=1.0))),
            ("hgb", make_pipeline(
                SimpleImputer(strategy="median"),
                HistGradientBoostingClassifier(random_state=SEED, max_iter=300,
                                               learning_rate=0.05))),
        ]:
            model.fit(Xtr, ytr)
            p = model.predict_proba(Xte)[:, 1]
            yhat = (p >= 0.5).astype(int)
            rows.append({
                "target": target, "split": split_col, "subset": subset_name,
                "feature_group": gname, "n_feats": len(feats), "model": mname,
                "n_train": len(ytr), "n_test": len(yte), "pos_rate_test": float(yte.mean()),
                "AUROC": roc_auc_score(yte, p),
                "AUPRC": average_precision_score(yte, p),
                "F1": f1_score(yte, yhat),
                "balACC": balanced_accuracy_score(yte, yhat),
                "ECE": ece(yte.astype(float), p),
            })
    return rows


def run_regression(target, split_col, subset_name, df):
    rows = []
    tr = df[df[split_col] == "train"]; te = df[df[split_col] == "test"]
    for gname, feats in groups.items():
        if not feats:
            continue
        Xtr, ytr = get_xy(tr, feats, target)
        Xte, yte = get_xy(te, feats, target)
        for mname, model in [
            ("ridge", make_pipeline(SimpleImputer(strategy="median"),
                                    StandardScaler(), Ridge(alpha=1.0))),
            ("hgb", make_pipeline(SimpleImputer(strategy="median"),
                                  HistGradientBoostingRegressor(random_state=SEED,
                                  max_iter=300, learning_rate=0.05))),
        ]:
            model.fit(Xtr, ytr)
            pred = model.predict(Xte)
            rmse = mean_squared_error(yte, pred) ** 0.5
            rows.append({
                "target": target, "split": split_col, "subset": subset_name,
                "feature_group": gname, "n_feats": len(feats), "model": mname,
                "pearson_r": float(stats.pearsonr(yte, pred)[0]),
                "spearman_rho": float(stats.spearmanr(yte, pred)[0]),
                "RMSE": float(rmse), "MAE": float(mean_absolute_error(yte, pred)),
            })
    return rows


def run_multiclass(split_col, subset_name, df):
    rows = []
    tr = df[df[split_col] == "train"]; te = df[df[split_col] == "test"]
    for gname, feats in groups.items():
        if not feats:
            continue
        Xtr, ytr = get_xy(tr, feats, "y_risk_class")
        Xte, yte = get_xy(te, feats, "y_risk_class")
        model = make_pipeline(SimpleImputer(strategy="median"),
                              HistGradientBoostingClassifier(random_state=SEED,
                              max_iter=300, learning_rate=0.05))
        model.fit(Xtr, ytr)
        yhat = model.predict(Xte)
        rows.append({
            "split": split_col, "subset": subset_name, "feature_group": gname,
            "model": "hgb",
            "macro_F1": f1_score(yte, yhat, average="macro"),
            "weighted_F1": f1_score(yte, yhat, average="weighted"),
            "balACC": balanced_accuracy_score(yte, yhat),
        })
    return rows


# ── Run all ─────────────────────────────────────────────────────────────────
subset_names = [s.strip() for s in args.subsets.split(",") if s.strip()]
split_cols = [s.strip() for s in args.splits.split(",") if s.strip()]
subsets = {}
for sname in subset_names:
    if sname == "all":
        subsets["all"] = master
    elif sname == "associated":
        subsets["associated"] = master[master["associated"]].copy()
    else:
        raise ValueError(f"Unknown subset: {sname}")

suffix = args.out_suffix
cls_rows, reg_rows, mc_rows = [], [], []
for sname, sdf in subsets.items():
    for split in split_cols:
        print(f"[{sname} | {split}] classification y_high_I2 ...")
        cls_rows += run_classification("y_high_I2", split, sname, sdf)
        print(f"[{sname} | {split}] regression I2 ...")
        reg_rows += run_regression("I2", split, sname, sdf)
        if not args.skip_riskclass:
            print(f"[{sname} | {split}] 3-class risk ...")
            mc_rows += run_multiclass(split, sname, sdf)
        gc.collect()

cls = pd.DataFrame(cls_rows); reg = pd.DataFrame(reg_rows); mc = pd.DataFrame(mc_rows)
cls.to_csv(out_dir / f"ablation_classification{suffix}.csv", index=False, float_format="%.4f")
reg.to_csv(out_dir / f"ablation_regression{suffix}.csv", index=False, float_format="%.4f")
mc.to_csv(out_dir / f"ablation_riskclass{suffix}.csv", index=False, float_format="%.4f")

# ── Summary ─────────────────────────────────────────────────────────────────
lines = []
def log(s): print(s); lines.append(s)

GROUP_ORDER = ["FST", "AF", "LD", "SEL", "AF_LD", "AF_SEL", "LD_SEL", "AF_LD_SEL"]
log("=" * 72)
log("MECHANISM ABLATION — y_high_I2 (I2>0.25), split_variant, HGB model")
log("=" * 72)
for sname in subsets:
    sub = cls[(cls.subset == sname) & (cls.split == "split_variant") & (cls.model == "hgb")]
    sub = sub.set_index("feature_group").reindex(GROUP_ORDER).dropna(subset=["AUROC"])
    log(f"\n[{sname}]  (test pos-rate ~{sub['pos_rate_test'].iloc[0]:.2f})")
    log(sub[["AUROC", "AUPRC", "F1", "balACC", "ECE"]].round(3).to_string())

log("\n" + "=" * 72)
log("I2 REGRESSION (Spearman rho), split_variant, HGB model")
log("=" * 72)
for sname in subsets:
    sub = reg[(reg.subset == sname) & (reg.split == "split_variant") & (reg.model == "hgb")]
    sub = sub.set_index("feature_group").reindex(GROUP_ORDER).dropna(subset=["pearson_r"])
    log(f"\n[{sname}]")
    log(sub[["pearson_r", "spearman_rho", "RMSE", "MAE"]].round(3).to_string())

(out_dir / f"ablation_summary{suffix}.txt").write_text("\n".join(lines), encoding="utf-8")
print(
    f"\nSaved ablation_classification{suffix}.csv / ablation_regression{suffix}.csv / "
    f"ablation_riskclass{suffix}.csv / ablation_summary{suffix}.txt"
)
