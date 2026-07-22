"""
Phase C1: Export PGS Catalog–compatible portability report card.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--pgs-id", default="PGS000018")
    p.add_argument(
        "--predictions",
        default=str(ROOT / "data/modeling/variant_portability_predictions.genomewide.parquet"),
    )
    p.add_argument(
        "--pgs-root",
        default=str(ROOT / "data/processed/pgs_grch38"),
    )
    p.add_argument(
        "--tiers",
        default=str(ROOT / "data/labels/finemap_tiers_genomewide.parquet"),
    )
    p.add_argument(
        "--shap",
        default=str(ROOT / "results/tables/shap_mechanism_attribution_genomewide.csv"),
        help="Optional per-variant dominant mechanism",
    )
    p.add_argument(
        "--out",
        default="",
        help="Default: results/tables/portability_report_card_{pgs}.tsv",
    )
    return p.parse_args()


def parse_vid(vid: str) -> tuple[str, str, str, str]:
    parts = str(vid).split(":")
    if len(parts) >= 4:
        return parts[0], parts[1], parts[2], parts[3]
    return parts[0] if parts else "", parts[1] if len(parts) > 1 else "", "", ""


def main() -> None:
    args = parse_args()
    pgs = args.pgs_id
    # Load harmonized weights
    weight_path = Path(args.pgs_root) / pgs / f"{pgs}.harmonized.tsv"
    if not weight_path.exists():
        # try any tsv
        cands = list((Path(args.pgs_root) / pgs).glob("*.tsv"))
        if not cands:
            raise SystemExit(f"Missing weights for {pgs}")
        weight_path = cands[0]
    w = pd.read_csv(weight_path, sep="\t", dtype=str)
    # normalize columns
    colmap = {c.lower(): c for c in w.columns}
    vid_col = colmap.get("variant_id") or colmap.get("id")
    ea_col = colmap.get("effect_allele") or colmap.get("allele1")
    wt_col = colmap.get("effect_weight") or colmap.get("weight") or colmap.get("beta")
    w = w.rename(columns={vid_col: "variant_id", ea_col: "effect_allele", wt_col: "effect_weight"})
    w["effect_weight"] = pd.to_numeric(w["effect_weight"], errors="coerce")

    # Predictions
    import pyarrow.dataset as ds

    dataset = ds.dataset(str(args.predictions), format="parquet")
    names = set(dataset.schema.names)
    risk_col = None
    for cand in [
        "predicted_risk",
        "y_prob",
        "pred_prob",
        "prob_high_I2",
        "portability_risk",
        "risk_prob",
        "prediction",
    ]:
        if cand in names:
            risk_col = cand
            break
    if risk_col is None:
        raise SystemExit(f"No risk column in {args.predictions}: {names}")
    cols = ["variant_id", risk_col]
    if "trait" in names:
        cols.append("trait")
    chunks = []
    for batch in dataset.scanner(columns=cols, batch_size=500_000).to_batches():
        chunks.append(batch.to_pandas())
    pred = pd.concat(chunks, ignore_index=True)
    risk = pred.groupby("variant_id", as_index=False)[risk_col].max()
    risk = risk.rename(columns={risk_col: "portability_risk_prob"})

    card = w.merge(risk, on="variant_id", how="left")
    # Fine-map tier
    tiers_path = Path(args.tiers)
    if tiers_path.exists():
        tiers = pd.read_parquet(tiers_path, columns=["variant_id", "finemap_tier"])
        tiers = tiers.drop_duplicates("variant_id")
        card = card.merge(tiers, on="variant_id", how="left")
    else:
        card["finemap_tier"] = np.nan

    shap_path = Path(args.shap)
    if shap_path.exists():
        sh = pd.read_csv(shap_path)
        if {"variant_id", "dominant_mechanism"}.issubset(sh.columns):
            card = card.merge(sh[["variant_id", "dominant_mechanism"]], on="variant_id", how="left")
    if "dominant_mechanism" not in card.columns:
        card["dominant_mechanism"] = "AF_LD_SEL"

    parsed = card["variant_id"].map(parse_vid)
    card["chr_name"] = parsed.map(lambda x: x[0])
    card["chr_position"] = parsed.map(lambda x: x[1])
    card["score_risk_percentile"] = card["portability_risk_prob"].rank(pct=True) * 100

    out_cols = [
        "variant_id",
        "chr_name",
        "chr_position",
        "effect_allele",
        "effect_weight",
        "portability_risk_prob",
        "dominant_mechanism",
        "finemap_tier",
        "score_risk_percentile",
    ]
    out = Path(args.out) if args.out else ROOT / f"results/tables/portability_report_card_{pgs}.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    card[out_cols].to_csv(out, sep="\t", index=False, float_format="%.6g")
    print(f"Saved {out} ({len(card):,} rows)")


if __name__ == "__main__":
    main()
