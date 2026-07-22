"""
Evaluation suite:
- bootstrap CIs for AUROC/AUPRC
- calibration curves
- label-permutation negative control
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.calibration import calibration_curve

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MASTER = ROOT / "data/modeling/master_variant_table.parquet"
DEFAULT_GROUPS = ROOT / "data/modeling/feature_groups.json"

SEED = 719
RNG = np.random.default_rng(SEED)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--master", default=str(DEFAULT_MASTER))
    p.add_argument("--groups", default=str(DEFAULT_GROUPS))
    p.add_argument("--subsets", default="all,associated", help="Comma list: all, associated")
    p.add_argument("--associated-only", action="store_true", help="Load only associated rows (saves RAM)")
    p.add_argument("--splits", default="split_variant,split_trait", help="Comma list of split columns")
    p.add_argument("--n-boot", type=int, default=100)
    p.add_argument("--out-suffix", default="", help="Suffix for outputs, e.g. .genomewide_partial6")
    return p.parse_args()


def ece(y, p, bins=10):
    edges = np.linspace(0, 1, bins + 1)
    out = 0.0
    for i in range(bins):
        m = (p >= edges[i]) & (p < edges[i + 1])
        if not m.any():
            continue
        out += m.mean() * abs(y[m].mean() - p[m].mean())
    return float(out)


def bootstrap_ci(y, p, n_boot=100):
    n = len(y)
    au, ap = [], []
    for _ in range(n_boot):
        idx = RNG.integers(0, n, n)
        yy = y[idx]
        pp = p[idx]
        if len(np.unique(yy)) < 2:
            continue
        au.append(roc_auc_score(yy, pp))
        ap.append(average_precision_score(yy, pp))
    if not au:
        return np.nan, np.nan, np.nan, np.nan
    return (
        float(np.percentile(au, 2.5)),
        float(np.percentile(au, 97.5)),
        float(np.percentile(ap, 2.5)),
        float(np.percentile(ap, 97.5)),
    )


def load_master(
    path: str,
    associated_only: bool,
    columns: list[str] | None = None,
) -> pd.DataFrame:
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


if __name__ == "__main__":
    args = parse_args()
    suffix = args.out_suffix
    OUT_TABLE = ROOT / f"results/tables/headline_metrics_ci{suffix}.csv"
    FIG_CAL = ROOT / f"results/figures/fig_calibration{suffix}.png"
    FIG_NEG = ROOT / f"results/figures/fig_negative_control{suffix}.png"
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    FIG_CAL.parent.mkdir(parents=True, exist_ok=True)

    groups = json.loads(Path(args.groups).read_text())
    feat_groups = {k: groups[k] for k in ["FST", "AF_LD", "AF_LD_SEL"] if k in groups}
    split_cols = [s.strip() for s in args.splits.split(",") if s.strip()]
    cols = ["associated", "y_high_I2"] + split_cols
    # source holdout needs source_group even when not listed as a split column name
    cols.append("source_group")
    for feats in feat_groups.values():
        cols.extend(feats)
    # unique preserve order
    seen: set[str] = set()
    cols = [c for c in cols if not (c in seen or seen.add(c))]
    d = load_master(args.master, args.associated_only, columns=cols)
    print(f"evaluate_suite: {len(d):,} rows from {args.master}")
    # drop missing optional columns after load
    if "source_group" not in d.columns:
        print("  note: source_group absent — skipping split_source")

    subset_names = [s.strip() for s in args.subsets.split(",") if s.strip()]
    subset_defs = []
    for sname in subset_names:
        if sname == "all":
            subset_defs.append(("all", d))
        elif sname == "associated":
            subset_defs.append(("associated", d[d["associated"]].copy()))
        else:
            raise ValueError(f"Unknown subset: {sname}")

    rows = []
    cal_curves = {}
    split_defs = [
        ("split_variant", lambda df: (df[df["split_variant"] == "train"], df[df["split_variant"] == "test"])),
        ("split_trait", lambda df: (df[df["split_trait"] == "train"], df[df["split_trait"] == "test"])),
    ]
    split_defs = [item for item in split_defs if item[0] in split_cols]
    if "source_group" in d.columns:
        split_defs.append(
            (
                "split_source",
                lambda df: (
                    df[df["source_group"] == "panukbb_only"],
                    df[df["source_group"] == "multisource"],
                ),
            )
        )

    for subset_name, subset_df in subset_defs:
        for split_col, splitter in split_defs:
            tr, te = splitter(subset_df)
            ytr = tr["y_high_I2"].astype(int).to_numpy()
            yte = te["y_high_I2"].astype(int).to_numpy()
            if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
                continue

            for gname, feats in feat_groups.items():
                print(f"[eval] {subset_name} {split_col} {gname}")
                Xtr = tr[feats].to_numpy(dtype=np.float32)
                Xte = te[feats].to_numpy(dtype=np.float32)
                model = Pipeline(
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
                    ],
                    memory=None,
                )
                model.fit(Xtr, ytr)
                p = model.predict_proba(Xte)[:, 1]
                au = roc_auc_score(yte, p)
                ap = average_precision_score(yte, p)
                au_lo, au_hi, ap_lo, ap_hi = bootstrap_ci(yte, p, n_boot=args.n_boot)
                key = (subset_name, split_col, gname)
                cal_curves[key] = calibration_curve(yte, p, n_bins=10, strategy="quantile")
                rows.append(
                    {
                        "subset": subset_name,
                        "split": split_col,
                        "feature_group": gname,
                        "AUROC": au,
                        "AUROC_lo": au_lo,
                        "AUROC_hi": au_hi,
                        "AUPRC": ap,
                        "AUPRC_lo": ap_lo,
                        "AUPRC_hi": ap_hi,
                        "ECE": ece(yte.astype(float), p),
                        "n_test": len(te),
                        "pos_rate": float(yte.mean()),
                    }
                )

                # Negative control: permute labels on training only.
                yperm = RNG.permutation(ytr)
                model.fit(Xtr, yperm)
                p_perm = model.predict_proba(Xte)[:, 1]
                au_p = roc_auc_score(yte, p_perm)
                ap_p = average_precision_score(yte, p_perm)
                au_lo, au_hi, ap_lo, ap_hi = bootstrap_ci(yte, p_perm, n_boot=args.n_boot)
                rows.append(
                    {
                        "subset": subset_name,
                        "split": split_col,
                        "feature_group": "PERMUTED",
                        "base_group": gname,
                        "AUROC": au_p,
                        "AUROC_lo": au_lo,
                        "AUROC_hi": au_hi,
                        "AUPRC": ap_p,
                        "AUPRC_lo": ap_lo,
                        "AUPRC_hi": ap_hi,
                        "ECE": ece(yte.astype(float), p_perm),
                        "n_test": len(te),
                        "pos_rate": float(yte.mean()),
                    }
                )

    out = pd.DataFrame(rows)
    out.to_csv(OUT_TABLE, index=False, float_format="%.4f")
    print(f"Saved {OUT_TABLE}")

    # Calibration figure (associated + split_variant).
    fig, ax = plt.subplots(figsize=(7, 5))
    for gname in ["FST", "AF_LD", "AF_LD_SEL"]:
        key = ("associated", "split_variant", gname)
        if key not in cal_curves:
            continue
        prob_true, prob_pred = cal_curves[key]
        ax.plot(prob_pred, prob_true, marker="o", label=gname)
    ax.plot([0, 1], [0, 1], "--", color="black", linewidth=1)
    ax.set_xlabel("Predicted probability")
    ax.set_ylabel("Observed frequency")
    ax.set_title("Calibration (associated, unseen-variant split)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_CAL, dpi=150)
    plt.close(fig)
    print(f"Saved {FIG_CAL}")

    # Negative control summary figure.
    sub = out[(out["subset"] == "associated") & (out["split"] == "split_variant")]
    real = sub[sub["feature_group"] != "PERMUTED"].set_index("feature_group")
    perm = sub[sub["feature_group"] == "PERMUTED"].groupby("base_group", as_index=True)["AUROC"].mean()
    xs = ["FST", "AF_LD", "AF_LD_SEL"]
    x = np.arange(len(xs))
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(x - 0.18, [real.loc[g, "AUROC"] if g in real.index else np.nan for g in xs], width=0.36, label="real")
    ax.bar(x + 0.18, [perm.get(g, np.nan) for g in xs], width=0.36, label="permuted")
    ax.axhline(0.5, linestyle="--", color="black", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(xs)
    ax.set_ylabel("AUROC")
    ax.set_title("Negative control (label permutation)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_NEG, dpi=150)
    plt.close(fig)
    print(f"Saved {FIG_NEG}")
