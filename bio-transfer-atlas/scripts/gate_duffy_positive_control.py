"""
Phase B4: Duffy-null (rs2814778 / ACKR1) neighborhood positive-control gate.

GRCh38 ACKR1 region ~chr1:159173968-159176290; Duffy-null rs2814778 ≈ 1:159204893.
Checks whether portability-risk predictions concentrate near this locus vs matched nulls.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DUFFY_CHR = "1"
DUFFY_POS = 159_204_893
WINDOW = 200_000


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--predictions",
        default=str(ROOT / "data/modeling/variant_portability_predictions.genomewide.parquet"),
    )
    p.add_argument(
        "--out",
        default=str(ROOT / "results/tables/duffy_positive_control_genomewide.csv"),
    )
    p.add_argument("--window-bp", type=int, default=WINDOW)
    p.add_argument("--n-null", type=int, default=50)
    p.add_argument("--seed", type=int, default=719)
    return p.parse_args()


def parse_vid(vid: str) -> tuple[str, int] | None:
    parts = str(vid).split(":")
    if len(parts) < 2:
        return None
    try:
        return parts[0], int(parts[1])
    except ValueError:
        return None


def main() -> None:
    args = parse_args()
    pred_path = Path(args.predictions)
    if not pred_path.exists():
        raise SystemExit(f"Missing predictions {pred_path}")

    # Stream only needed cols
    import pyarrow.dataset as ds

    dataset = ds.dataset(str(pred_path), format="parquet")
    names = set(dataset.schema.names)
    risk_col = None
    for cand in [
        "predicted_risk",
        "portability_risk",
        "risk_prob",
        "y_prob",
        "pred_prob",
        "prob_high_I2",
        "risk",
        "prediction",
    ]:
        if cand in names:
            risk_col = cand
            break
    cols = ["variant_id"]
    if "trait" in names:
        cols.append("trait")
    if risk_col:
        cols.append(risk_col)
    else:
        # fall back: any float col
        for n in names:
            if n not in ("variant_id", "trait"):
                risk_col = n
                cols.append(n)
                break
    print(f"Using risk column: {risk_col}", flush=True)
    chunks = []
    for batch in dataset.scanner(columns=cols, batch_size=500_000).to_batches():
        chunks.append(batch.to_pandas())
    pred = pd.concat(chunks, ignore_index=True)
    if risk_col not in pred.columns:
        raise SystemExit(f"Could not find risk column in {pred_path}")

    # Aggregate max risk per variant across traits
    g = pred.groupby("variant_id", as_index=False)[risk_col].max()
    coords = g["variant_id"].map(parse_vid)
    g["chrom"] = coords.map(lambda x: x[0] if x else None)
    g["pos"] = coords.map(lambda x: x[1] if x else np.nan)
    g = g.dropna(subset=["pos"])

    duffy = g[
        (g["chrom"] == DUFFY_CHR)
        & (g["pos"] >= DUFFY_POS - args.window_bp)
        & (g["pos"] <= DUFFY_POS + args.window_bp)
    ]
    # Matched null: same chrom length windows elsewhere on chr1
    chr1 = g[g["chrom"] == DUFFY_CHR].copy()
    rng = np.random.default_rng(args.seed)
    null_means = []
    positions = chr1["pos"].to_numpy()
    lo, hi = int(positions.min()), int(positions.max())
    for _ in range(args.n_null):
        center = int(rng.integers(lo + args.window_bp, hi - args.window_bp))
        if abs(center - DUFFY_POS) < 2 * args.window_bp:
            continue
        win = chr1[(chr1["pos"] >= center - args.window_bp) & (chr1["pos"] <= center + args.window_bp)]
        if len(win):
            null_means.append(float(win[risk_col].mean()))

    duffy_mean = float(duffy[risk_col].mean()) if len(duffy) else float("nan")
    null_mean = float(np.mean(null_means)) if null_means else float("nan")
    null_sd = float(np.std(null_means)) if null_means else float("nan")
    z = (duffy_mean - null_mean) / null_sd if null_sd and null_sd > 0 else float("nan")
    # Genome-wide percentile of duffy neighborhood mean
    all_means_proxy = float(g[risk_col].mean())
    pct = float((g[risk_col] < duffy_mean).mean()) if len(duffy) else float("nan")

    row = {
        "locus": "ACKR1_Duffy_null_rs2814778",
        "chrom": DUFFY_CHR,
        "pos": DUFFY_POS,
        "window_bp": args.window_bp,
        "n_variants_window": len(duffy),
        "mean_risk_window": duffy_mean,
        "mean_risk_null_windows": null_mean,
        "sd_risk_null_windows": null_sd,
        "z_vs_null": z,
        "frac_genome_below_duffy_mean": pct,
        "genome_mean_risk": all_means_proxy,
        "passes_positive_control": bool(z > 1.0) if pd.notna(z) else False,
    }
    out = pd.DataFrame([row])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, float_format="%.6g")
    print(out.to_string(index=False))
    print(f"Saved {args.out}")
    if not row["passes_positive_control"]:
        print(
            "NOTE: Duffy neighborhood did not exceed null z>1 — "
            "expected if WBC trait not yet in training labels; gate records result honestly."
        )


if __name__ == "__main__":
    main()
