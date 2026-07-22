#!/usr/bin/env python3
"""
Compute autosomal (or available-chrom) Z-correlation rg from existing Pan-UKB chrom parquets.
Faster companion while full Popcorn downloads run.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/features/baselines"
TRAITS = ["T2D", "CAD", "BMI", "LDL"]
PAIRS = [("EUR", "AFR"), ("EUR", "EAS")]


def z_from(df, anc):
    b, s = f"beta_{anc}", f"se_{anc}"
    if b in df.columns and s in df.columns:
        se = pd.to_numeric(df[s], errors="coerce").replace(0, np.nan)
        return pd.to_numeric(df[b], errors="coerce") / se
    return pd.to_numeric(df[b], errors="coerce")


def main():
    avail = []
    for c in range(1, 23):
        d = ROOT / "data/raw/panukbb" / f"chr{c}"
        if any((d / f"{t}.chr{c}.parquet").exists() for t in TRAITS):
            avail.append(str(c))
    print("available chroms", avail, flush=True)
    rows = []
    for trait in TRAITS:
        frames = []
        for c in avail:
            p = ROOT / "data/raw/panukbb" / f"chr{c}" / f"{trait}.chr{c}.parquet"
            if p.exists():
                frames.append(pd.read_parquet(p))
        if not frames:
            continue
        df = pd.concat(frames, ignore_index=True)
        print(f"{trait}: {len(df):,} rows", flush=True)
        for a1, a2 in PAIRS:
            z1 = z_from(df, a1).to_numpy(float)
            z2 = z_from(df, a2).to_numpy(float)
            m = np.isfinite(z1) & np.isfinite(z2)
            n = int(m.sum())
            rg = float(np.corrcoef(z1[m], z2[m])[0, 1]) if n > 200 else float("nan")
            se = float(1.0 / np.sqrt(max(n - 3, 1)))
            method = (
                "panukbb_z_concordance_autosome"
                if len(avail) >= 20
                else f"panukbb_z_concordance_chr{','.join(avail)}"
            )
            rows.append(
                {
                    "trait": trait,
                    "anc1": a1,
                    "anc2": a2,
                    "z_concordance": rg,
                    "rg": rg,  # legacy alias
                    "se": se,
                    "n": n,
                    "method": method,
                    "estimand": "cross_ancestry_z_score_concordance",
                }
            )
            print(f"  {a1}-{a2} z_concordance={rg:.4f} n={n}", flush=True)
    summary = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    (ROOT / "results/tables").mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT / "popcorn_rg_summary.csv", index=False, float_format="%.6g")
    summary.to_csv(ROOT / "results/tables/popcorn_rg_summary.csv", index=False, float_format="%.6g")
    summary.to_csv(ROOT / "results/tables/z_concordance_by_trait_pair.csv", index=False, float_format="%.6g")
    summary.to_csv(ROOT / "results/tables/ldsc_rg_companion.csv", index=False, float_format="%.6g")
    feat = []
    for trait in TRAITS:
        rec = {"trait": trait}
        for a1, a2 in PAIRS:
            hit = summary[(summary.trait == trait) & (summary.anc1 == a1) & (summary.anc2 == a2)]
            zc = float(hit["z_concordance"].iloc[0]) if len(hit) else np.nan
            se_v = float(hit["se"].iloc[0]) if len(hit) else np.nan
            rec[f"z_concordance_{a1}_{a2}"] = zc
            rec[f"z_concordance_{a1}_{a2}_se"] = se_v
            rec[f"rg_{a1}_{a2}"] = zc  # legacy
            rec[f"rg_{a1}_{a2}_se"] = se_v
        feat.append(rec)
    pd.DataFrame(feat).to_parquet(OUT / "rg_real_by_trait.parquet", index=False)
    pd.DataFrame(feat).to_parquet(OUT / "z_concordance_by_trait.parquet", index=False)
    meta = {
        "chroms": avail,
        "method": summary["method"].iloc[0] if len(summary) else "none",
        "scope": "available_chroms",
        "estimand": "cross_ancestry_z_score_concordance",
        "note": "Pearson Z-score concordance; not Popcorn genetic-effect rg unless popcorn fit succeeds",
    }
    (OUT / "rg_real_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    (ROOT / "results/tables/rg_real_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("Z_CONCORDANCE_AVAILABLE_DONE")


if __name__ == "__main__":
    main()
