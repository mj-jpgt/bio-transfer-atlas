"""
Phase B3: VEP–ancestry interaction features.

Joins AlphaMissense if present; otherwise uses LOEUF/pLI × AF divergence as proxy.
Writes data/features/selection/vep_af_interaction_features.parquet and a stratum eval table.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
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
        "--alphamissense",
        default=str(ROOT / "data/annotations/alphamissense_grch38.parquet"),
        help="Optional: variant_id, alphamissense_score",
    )
    p.add_argument(
        "--groups",
        default=str(ROOT / "data/modeling/feature_groups_genomewide_genomewide.json"),
    )
    p.add_argument(
        "--out-features",
        default=str(ROOT / "data/features/selection/vep_af_interaction_features.parquet"),
    )
    p.add_argument(
        "--out-eval",
        default=str(ROOT / "results/tables/vep_af_interaction_eval.csv"),
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


def main() -> None:
    args = parse_args()
    cols = [
        "variant_id",
        "y_high_I2",
        "split_variant",
        "associated",
        "AF_max_diff",
        "AF_diff_AFR_EUR",
        "AF_EUR",
        "LOEUF",
        "pLI",
        "mis_z",
    ]
    print("Streaming sample ...", flush=True)
    # Need AF_LD_SEL cols for eval too
    groups_preview = json.loads(Path(args.groups).read_text(encoding="utf-8"))
    cols = list(
        dict.fromkeys(
            cols
            + groups_preview.get("AF_LD_SEL", [])
        )
    )
    sample_path = ROOT / "data/modeling/_tmp_vep_eval_sample.parquet"
    stream_sample_associated(
        Path(args.master), cols, sample_path, max_rows=args.max_sample, seed=SEED
    )
    df = pd.read_parquet(sample_path)
    print(f"Sampled {len(df):,}", flush=True)
    # Unique variants for feature table
    feats = (
        df.groupby("variant_id", as_index=False)
        .agg(
            AF_max_diff=("AF_max_diff", "first"),
            AF_diff_AFR_EUR=("AF_diff_AFR_EUR", "first"),
            AF_EUR=("AF_EUR", "first"),
            LOEUF=("LOEUF", "first"),
            pLI=("pLI", "first"),
            mis_z=("mis_z", "first"),
        )
    )
    am_path = Path(args.alphamissense)
    if am_path.exists():
        am = pd.read_parquet(am_path)
        feats = feats.merge(am, on="variant_id", how="left")
        score_col = "alphamissense_score"
        feats["vep_source"] = "alphamissense"
    else:
        # Constraint severity proxy (higher = more constrained / severe)
        feats["alphamissense_score"] = (
            (1.0 - feats["LOEUF"].clip(0, 1).fillna(0.5)) * 0.5
            + feats["pLI"].fillna(0).clip(0, 1) * 0.5
        )
        score_col = "alphamissense_score"
        feats["vep_source"] = "loeuf_pli_proxy"

    feats["vep_x_afmax"] = feats[score_col] * feats["AF_max_diff"].fillna(0)
    feats["vep_x_afr_eur"] = feats[score_col] * feats["AF_diff_AFR_EUR"].fillna(0)
    # Rare-het-like stratum: low EUR AF
    feats["rare_eur"] = (feats["AF_EUR"].fillna(0.5) < 0.01).astype(float)

    Path(args.out_features).parent.mkdir(parents=True, exist_ok=True)
    feats.to_parquet(args.out_features, index=False)
    print(f"Saved {args.out_features} ({len(feats):,}) source={feats['vep_source'].iloc[0]}")

    # Eval: AF_LD_SEL vs + interaction on associated rows
    groups = json.loads(Path(args.groups).read_text(encoding="utf-8"))
    base_feats = [c for c in groups["AF_LD_SEL"] if c in df.columns]
    df = df.merge(
        feats[["variant_id", "vep_x_afmax", "vep_x_afr_eur", "rare_eur", score_col]],
        on="variant_id",
        how="left",
    )

    def run(feat_list: list[str], name: str, mask: pd.Series | None = None) -> dict:
        import gc

        sub = df if mask is None else df[mask]
        tr = cap_train(sub[sub["split_variant"] == "train"], args.max_train)
        te = sub[sub["split_variant"] == "test"]
        use = [f for f in feat_list if f in sub.columns]
        ytr = tr["y_high_I2"].astype(int).to_numpy()
        yte = te["y_high_I2"].astype(int).to_numpy()
        if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
            return {"setting": name, "AUROC": np.nan, "n_test": len(te)}
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
        Xtr = tr[use].to_numpy(dtype=np.float32, copy=True)
        Xte = te[use].to_numpy(dtype=np.float32, copy=True)
        model.fit(Xtr, ytr)
        del Xtr
        gc.collect()
        p = model.predict_proba(Xte)[:, 1]
        del Xte, model
        gc.collect()
        return {"setting": name, "AUROC": float(roc_auc_score(yte, p)), "n_test": len(te)}

    rows = [
        run(base_feats, "AF_LD_SEL"),
        run(base_feats + ["vep_x_afmax", "vep_x_afr_eur", score_col], "AF_LD_SEL+VEP_AF"),
        run(
            base_feats + ["vep_x_afmax", "vep_x_afr_eur", score_col],
            "AF_LD_SEL+VEP_AF_rare_eur",
            mask=df["rare_eur"] == 1,
        ),
        run(
            base_feats + ["vep_x_afmax", "vep_x_afr_eur", score_col],
            "AF_LD_SEL+VEP_AF_common_eur",
            mask=df["rare_eur"] == 0,
        ),
    ]
    out = pd.DataFrame(rows)
    Path(args.out_eval).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_eval, index=False, float_format="%.4f")
    print(out.to_string(index=False))
    print(f"Saved {args.out_eval}")

    # Extend feature groups
    g = groups
    g["VEP_AF"] = [score_col, "vep_x_afmax", "vep_x_afr_eur", "rare_eur"]
    Path(args.groups).write_text(json.dumps(g, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
