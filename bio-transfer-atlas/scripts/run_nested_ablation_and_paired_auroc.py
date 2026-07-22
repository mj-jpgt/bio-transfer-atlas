#!/usr/bin/env python3
"""Nested AF/LD/SEL ablations + paired LD-block bootstrap ΔAUROC/AUPRC with verdicts."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
SEED = 719


def auroc_safe(y, p):
    y = np.asarray(y).astype(int)
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, p))


def auprc_safe(y, p):
    y = np.asarray(y).astype(int)
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(average_precision_score(y, p))


def verdict_from_ci(lo: float, hi: float) -> str:
    if not np.isfinite(lo) or not np.isfinite(hi):
        return "inconclusive"
    if lo > 0:
        return "higher"
    if hi < 0:
        return "lower"
    return "no_clear_difference"


def nest_feature_lists(groups: dict, columns: set[str]) -> dict[str, list[str]]:
    """Prefer explicit JSON keys; never put FST into AF."""

    def clean(cols: list[str]) -> list[str]:
        return [c for c in cols if c in columns and not (c.startswith("FST_") or c == "FST_like")]

    af = clean(groups.get("AF", []))
    ld = clean(groups.get("LD", []))
    sel = clean(groups.get("SEL", []))
    # Fallback heuristics without FST contamination
    if not af or not ld:
        all_af = clean(groups.get("AF_LD_SEL", []))
        if not af:
            af = [f for f in all_af if f.startswith("AF_") or f == "AIS" or "freq" in f.lower()]
        if not ld:
            ld = [f for f in all_af if f.startswith("LD_") or f == "LD_entropy"]
        if not sel:
            sel = [
                f
                for f in all_af
                if f not in af
                and f not in ld
                and not f.startswith("FST_")
                and f != "FST_like"
            ]
    af_ld_sel = clean(groups.get("AF_LD_SEL", [])) or sorted(set(af + ld + sel))
    return {
        "AF": af,
        "LD": ld,
        "SEL": sel,
        "AF_LD": sorted(set(af + ld)),
        "AF_SEL": sorted(set(af + sel)),
        "LD_SEL": sorted(set(ld + sel)),
        "AF_LD_SEL": af_ld_sel,
    }


def main() -> None:
    sample = ROOT / "data/modeling/_tmp_ldblock_associated_sample.parquet"
    groups = json.loads(
        (ROOT / "data/modeling/feature_groups_genomewide_genomewide.json").read_text(encoding="utf-8")
    )
    df = pd.read_parquet(sample)
    blocks = ROOT / "data/modeling/ld_block_assignments_genomewide.parquet"
    if "split_ld_block" not in df.columns and blocks.exists():
        ld = pd.read_parquet(blocks, columns=["variant_id", "split_ld_block", "ld_block"])
        df = df.merge(ld.drop_duplicates("variant_id"), on="variant_id", how="inner")
    if "ld_block" not in df.columns and blocks.exists():
        ld = pd.read_parquet(blocks, columns=["variant_id", "ld_block"])
        df = df.merge(ld.drop_duplicates("variant_id"), on="variant_id", how="left")

    nest = nest_feature_lists(groups, set(df.columns))
    pop = [f for f in groups.get("POP_DISTANCE", []) if f in df.columns]
    fst = [f for f in groups.get("FST", ["FST_like"]) if f in df.columns]

    variant_feature_columns = sorted(
        set(sum((nest[k] for k in nest), [])) | set(pop) | set(fst)
    )
    assert not any(
        "z_concordance" in c or c.startswith("rg_") for c in variant_feature_columns
    ), "trait-constant concordance leaked into variant feature contest"
    assert not any(c.startswith("FST_") or c == "FST_like" for c in nest["AF"]), (
        "FST leaked into AF nest"
    )

    split = "split_ld_block" if "split_ld_block" in df.columns else None
    if split:
        tr_mask = df[split] == "train"
        te_mask = df[split] == "test"
    else:
        rng = np.random.default_rng(SEED)
        te_mask = pd.Series(rng.random(len(df)) < 0.2, index=df.index)
        tr_mask = ~te_mask

    ytr = df.loc[tr_mask, "y_high_I2"].astype(int)
    yte = df.loc[te_mask, "y_high_I2"].astype(int)
    rows = []
    preds = {}
    for name, feats in {**nest, "POP_DISTANCE": pop, "FST": fst}.items():
        feats = [f for f in feats if f in df.columns]
        if len(feats) < 1:
            continue
        imp = SimpleImputer(strategy="median")
        Xtr = imp.fit_transform(df.loc[tr_mask, feats])
        Xte = imp.transform(df.loc[te_mask, feats])
        clf = HistGradientBoostingClassifier(
            max_depth=6, learning_rate=0.05, max_iter=200, random_state=SEED
        )
        clf.fit(Xtr, ytr)
        p = clf.predict_proba(Xte)[:, 1]
        preds[name] = p
        rows.append(
            {
                "feature_group": name,
                "split": "split_ld_block",
                "AUROC": auroc_safe(yte, p),
                "AUPRC": auprc_safe(yte, p),
                "n_feats": len(feats),
                "n_test": int(te_mask.sum()),
            }
        )
        print(name, rows[-1]["AUROC"], rows[-1]["AUPRC"], flush=True)

    nest_path = ROOT / "results/tables/ablation_nested_af_ld_sel.csv"
    pd.DataFrame(rows).to_csv(nest_path, index=False, float_format="%.4f")

    te = df.loc[te_mask].copy()
    te["_idx"] = np.arange(len(te))
    block_col = "ld_block" if "ld_block" in te.columns else None
    rng = np.random.default_rng(SEED)
    pair_keys = [
        ("AF_LD_SEL", "POP_DISTANCE"),
        ("AF_LD_SEL", "FST"),
        ("AF_LD_SEL", "AF"),
        ("AF_LD", "AF"),
        ("AF_LD_SEL", "AF_LD"),
    ]
    deltas_auroc = {k: [] for k in pair_keys}
    deltas_auprc = {k: [] for k in pair_keys}
    if block_col and "AF_LD_SEL" in preds:
        blocks_u = te[block_col].unique()
        yte_np = yte.to_numpy()
        for _ in range(200):
            samp_blocks = rng.choice(blocks_u, size=len(blocks_u), replace=True)
            idx = te[te[block_col].isin(samp_blocks)]["_idx"].to_numpy()
            if len(idx) < 100:
                continue
            yb = yte_np[idx]
            for a, b in pair_keys:
                if a not in preds or b not in preds:
                    continue
                da = auroc_safe(yb, preds[a][idx])
                db = auroc_safe(yb, preds[b][idx])
                if np.isfinite(da) and np.isfinite(db):
                    deltas_auroc[(a, b)].append(da - db)
                pa = auprc_safe(yb, preds[a][idx])
                pb = auprc_safe(yb, preds[b][idx])
                if np.isfinite(pa) and np.isfinite(pb):
                    deltas_auprc[(a, b)].append(pa - pb)

    auroc_map = {r["feature_group"]: r["AUROC"] for r in rows}
    auprc_map = {r["feature_group"]: r["AUPRC"] for r in rows}

    def summarize(metric: str, point_map: dict, boots: dict) -> list[dict]:
        out = []
        for (a, b), vals in boots.items():
            if a not in point_map or b not in point_map:
                continue
            point = point_map[a] - point_map[b]
            if not vals:
                out.append(
                    {
                        "metric": metric,
                        "model_a": a,
                        "model_b": b,
                        "delta": point,
                        "delta_lo": np.nan,
                        "delta_hi": np.nan,
                        "p_gt_0": np.nan,
                        "n_boot": 0,
                        "ci_excludes_zero": False,
                        "verdict": "inconclusive",
                        "status": "point_only",
                    }
                )
                continue
            lo, hi = float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
            p_gt = float(np.mean(np.asarray(vals) > 0))
            out.append(
                {
                    "metric": metric,
                    "model_a": a,
                    "model_b": b,
                    "delta": point,
                    "delta_lo": lo,
                    "delta_hi": hi,
                    "p_gt_0": p_gt,
                    "n_boot": len(vals),
                    "ci_excludes_zero": bool(lo > 0 or hi < 0),
                    "verdict": verdict_from_ci(lo, hi),
                    "status": "ldblock_bootstrap",
                }
            )
        return out

    pair_rows = summarize("AUROC", auroc_map, deltas_auroc)
    # legacy columns for gate compatibility
    for r in pair_rows:
        r["delta_auroc"] = r["delta"]
    auprc_rows = summarize("AUPRC", auprc_map, deltas_auprc)

    pair_path = ROOT / "results/tables/auroc_paired_delta_ldblock.csv"
    pd.DataFrame(pair_rows).to_csv(pair_path, index=False, float_format="%.4f")
    auprc_path = ROOT / "results/tables/auprc_paired_delta_ldblock.csv"
    pd.DataFrame(auprc_rows).to_csv(auprc_path, index=False, float_format="%.4f")
    print(f"Saved {nest_path}, {pair_path}, {auprc_path}")


if __name__ == "__main__":
    main()
