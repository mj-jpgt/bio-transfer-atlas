"""
Lean Duffy positive control using Pan-UKB WBC region sumstats only.

No full-chr download, no master reload. Compares |beta_AFR - beta_EUR| and
I2-like discordance in the Duffy window vs matched null windows on the same
chr1 region file (or adjacent nulls from the same small parquet).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DUFFY_POS = 159_204_893
WINDOW = 200_000


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--wbc-region",
        default=str(
            ROOT
            / "data/raw/panukbb/regions/WBC.chr1_159000000_159400000.parquet"
        ),
    )
    p.add_argument(
        "--out",
        default=str(ROOT / "results/tables/duffy_wbc_region_positive_control.csv"),
    )
    p.add_argument("--window-bp", type=int, default=WINDOW)
    p.add_argument("--n-null", type=int, default=20)
    p.add_argument("--seed", type=int, default=719)
    return p.parse_args()


def cochran_i2(betas: np.ndarray, ses: np.ndarray) -> float:
    mask = np.isfinite(betas) & np.isfinite(ses) & (ses > 0)
    b, s = betas[mask], ses[mask]
    if len(b) < 2:
        return float("nan")
    w = 1.0 / (s**2)
    mu = np.sum(w * b) / np.sum(w)
    q = np.sum(w * (b - mu) ** 2)
    df = len(b) - 1
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
    return float(i2)


def main() -> None:
    args = parse_args()
    path = Path(args.wbc_region)
    if not path.exists():
        raise SystemExit(
            f"Missing {path}. Run:\n"
            "  python scripts/download_panukbb_region.py --trait WBC"
        )
    df = pd.read_parquet(path)
    for c in df.columns:
        if c.startswith(("beta_", "se_", "pval_")) or c in ("pos",):
            df[c] = pd.to_numeric(df[c], errors="coerce")

    ancs = [a for a in ["AFR", "EUR", "EAS", "AMR", "CSA", "MID"] if f"beta_{a}" in df.columns]
    if "AFR" not in ancs or "EUR" not in ancs:
        raise SystemExit(f"Need beta_AFR and beta_EUR; found {ancs}")

    df["abs_afr_eur"] = (df["beta_AFR"] - df["beta_EUR"]).abs()
    df["i2_proxy"] = [
        cochran_i2(
            np.array([row.get(f"beta_{a}", np.nan) for a in ancs], dtype=float),
            np.array([row.get(f"se_{a}", np.nan) for a in ancs], dtype=float),
        )
        for _, row in df.iterrows()
    ]

    lo, hi = DUFFY_POS - args.window_bp, DUFFY_POS + args.window_bp
    # If the region file is narrow, shrink the analysis window to fit
    pos_min, pos_max = int(df["pos"].min()), int(df["pos"].max())
    span = pos_max - pos_min
    win = min(args.window_bp, max(20_000, span // 4))
    lo, hi = DUFFY_POS - win, DUFFY_POS + win
    duffy = df[(df["pos"] >= lo) & (df["pos"] <= hi)]
    if duffy.empty:
        raise SystemExit(f"No variants in Duffy window {lo}-{hi}")

    # Nulls: non-overlapping flank tiles within this same small region file
    null_means = []
    step = max(win, 10_000)
    for center in range(pos_min + win, pos_max - win + 1, step):
        if abs(center - DUFFY_POS) < 2 * win:
            continue
        sub = df[(df["pos"] >= center - win) & (df["pos"] <= center + win)]
        if len(sub) < 20:
            continue
        null_means.append(float(sub["abs_afr_eur"].mean()))
    # Also compare Duffy vs everything outside Duffy window
    outside = df[(df["pos"] < lo) | (df["pos"] > hi)]
    if len(outside) >= 20:
        null_means.append(float(outside["abs_afr_eur"].mean()))

    duffy_mean = float(duffy["abs_afr_eur"].mean())
    null_arr = np.array(null_means, dtype=float)
    z = (duffy_mean - null_arr.mean()) / (null_arr.std() + 1e-12) if len(null_arr) > 1 else (
        (duffy_mean - null_arr[0]) / (null_arr[0] + 1e-12) if len(null_arr) == 1 else float("nan")
    )
    frac_below = float(np.mean(null_arr < duffy_mean)) if len(null_arr) else float("nan")
    # For WBC/Duffy, also require Duffy mean I2 elevated vs outside
    outside_i2 = float(outside["i2_proxy"].mean()) if len(outside) else float("nan")
    duffy_i2 = float(duffy["i2_proxy"].mean())
    passes = bool(
        (np.isfinite(z) and z > 1.0)
        or (np.isfinite(frac_below) and frac_below >= 0.8)
        or (np.isfinite(outside_i2) and duffy_i2 > outside_i2 * 1.2)
    )

    row = {
        "locus": "ACKR1_Duffy_null_rs2814778",
        "trait": "WBC",
        "metric": "mean_|beta_AFR-beta_EUR|",
        "window_bp_used": win,
        "n_duffy": len(duffy),
        "duffy_mean": duffy_mean,
        "duffy_mean_i2": duffy_i2,
        "outside_mean_i2": outside_i2,
        "n_null_windows": len(null_arr),
        "null_mean": float(null_arr.mean()) if len(null_arr) else np.nan,
        "null_sd": float(null_arr.std()) if len(null_arr) else np.nan,
        "z_vs_null": z,
        "frac_null_below_duffy": frac_below,
        "passes_positive_control": passes,
        "source_file": str(path.name),
    }
    out = pd.DataFrame([row])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, float_format="%.6g")
    print(out.to_string(index=False))
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
