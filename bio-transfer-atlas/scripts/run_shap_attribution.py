"""
Phase D: SHAP mechanism attribution on HGB (primary interpretability).
Optional lightweight LD-graph GNN smoke test if torch_geometric is installed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

ROOT = Path(__file__).resolve().parents[1]
SEED = 719


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
    p.add_argument("--max-rows", type=int, default=80_000, help="Subsample for SHAP runtime")
    p.add_argument(
        "--out",
        default=str(ROOT / "results/tables/shap_mechanism_attribution_genomewide.csv"),
    )
    p.add_argument("--try-gnn", action="store_true")
    return p.parse_args()


def load_sample(master: Path, cols: list[str], max_rows: int) -> pd.DataFrame:
    import pyarrow.dataset as ds

    dataset = ds.dataset(str(master), format="parquet")
    use = [c for c in cols if c in set(dataset.schema.names)]
    filt = ds.field("associated") == True  # noqa: E712
    chunks = []
    n = 0
    for batch in dataset.scanner(columns=use, filter=filt, batch_size=200_000).to_batches():
        chunks.append(batch.to_pandas())
        n += len(chunks[-1])
        if n >= max_rows * 3:
            break
    df = pd.concat(chunks, ignore_index=True)
    if len(df) > max_rows:
        df = df.sample(max_rows, random_state=SEED)
    return df


def mechanism_from_feats(feat_names: list[str], abs_shap: np.ndarray) -> str:
    groups = {"AF": 0.0, "LD": 0.0, "SEL": 0.0, "OTHER": 0.0}
    for name, val in zip(feat_names, abs_shap):
        if name.startswith("AF_") or name in ("FST_like", "AIS"):
            groups["AF"] += val
        elif name.startswith("LD_"):
            groups["LD"] += val
        elif name.startswith(("PBS", "FST_", "LOEUF", "pLI", "mis_z", "DAF")):
            groups["SEL"] += val
        else:
            groups["OTHER"] += val
    return max(groups, key=groups.get)


def main() -> None:
    args = parse_args()
    groups = json.loads(Path(args.groups).read_text(encoding="utf-8"))
    feats = groups["AF_LD_SEL"]
    cols = ["variant_id", "y_high_I2", "split_variant"] + feats
    print("Loading sample ...", flush=True)
    df = load_sample(Path(args.master), cols, args.max_rows)
    use = [f for f in feats if f in df.columns]
    tr = df[df["split_variant"] == "train"]
    te = df[df["split_variant"] == "test"]
    if te.empty:
        te = df.sample(min(5000, len(df)), random_state=SEED)
        tr = df.drop(te.index)

    pipe = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            (
                "hgb",
                HistGradientBoostingClassifier(
                    random_state=SEED, max_iter=150, learning_rate=0.05
                ),
            ),
        ]
    )
    Xtr = tr[use].to_numpy(dtype=np.float32)
    ytr = tr["y_high_I2"].astype(int).to_numpy()
    pipe.fit(Xtr, ytr)
    Xte = te[use].to_numpy(dtype=np.float32)

    # Prefer shap; fallback to permutation importance on a small set
    rows = []
    try:
        import shap

        # TreeExplainer needs the HGB step; use imputed matrix
        Xte_imp = pipe.named_steps["imputer"].transform(Xte)
        explainer = shap.Explainer(pipe.named_steps["hgb"].predict_proba, Xte_imp[:200])
        sv = explainer(Xte_imp[:2000])
        # class 1
        vals = sv.values
        if vals.ndim == 3:
            vals = vals[:, :, 1]
        abs_mean = np.abs(vals).mean(axis=0)
        for i, vid in enumerate(te["variant_id"].iloc[:2000]):
            dom = mechanism_from_feats(use, np.abs(vals[i]))
            rows.append(
                {
                    "variant_id": vid,
                    "dominant_mechanism": dom,
                    "shap_af": float(
                        sum(np.abs(vals[i][j]) for j, n in enumerate(use) if n.startswith("AF_") or n in ("FST_like", "AIS"))
                    ),
                    "shap_ld": float(
                        sum(np.abs(vals[i][j]) for j, n in enumerate(use) if n.startswith("LD_"))
                    ),
                    "shap_sel": float(
                        sum(
                            np.abs(vals[i][j])
                            for j, n in enumerate(use)
                            if n.startswith(("PBS", "FST_", "LOEUF", "pLI", "mis_z", "DAF"))
                        )
                    ),
                }
            )
        # global summary
        global_dom = mechanism_from_feats(use, abs_mean)
        print(f"SHAP global dominant mechanism: {global_dom}")
        method = "shap"
    except Exception as e:
        print(f"SHAP unavailable ({e}); using feature-group permutation proxy", flush=True)
        from sklearn.inspection import permutation_importance

        Xte_imp = pipe.named_steps["imputer"].transform(Xte[:3000])
        yte = te["y_high_I2"].astype(int).to_numpy()[:3000]
        r = permutation_importance(
            pipe.named_steps["hgb"], Xte_imp, yte, n_repeats=3, random_state=SEED, scoring="roc_auc"
        )
        abs_imp = np.maximum(r.importances_mean, 0)
        for i, vid in enumerate(te["variant_id"].iloc[:3000]):
            # assign same global mechanism (variant-level approx)
            rows.append(
                {
                    "variant_id": vid,
                    "dominant_mechanism": mechanism_from_feats(use, abs_imp),
                }
            )
        method = "permutation_importance"

    out = pd.DataFrame(rows)
    out["attribution_method"] = method
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, float_format="%.6g")
    print(out["dominant_mechanism"].value_counts())
    print(f"Saved {args.out}")

    if args.try_gnn:
        try:
            import torch
            from torch import nn

            print(f"torch {torch.__version__} cuda={torch.cuda.is_available()}")
            # Minimal MLP baseline standing in for GNN when torch_geometric absent
            try:
                import torch_geometric  # noqa: F401

                print("torch_geometric available — full GAT deferred to dedicated trainer")
            except ImportError:
                print("torch_geometric not installed; recorded MLP-ready torch env for Phase D")
            (ROOT / "results/tables/gnn_phase_d_status.txt").write_text(
                f"torch_ok cuda={torch.cuda.is_available()} geometric="
                + str(_has_geometric())
                + "\nHGB+SHAP is primary attribution; GAT optional follow-on.\n",
                encoding="utf-8",
            )
        except Exception as e:
            print(f"GNN smoke skipped: {e}")


def _has_geometric() -> bool:
    try:
        import torch_geometric  # noqa: F401

        return True
    except ImportError:
        return False


if __name__ == "__main__":
    main()
