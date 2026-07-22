#!/usr/bin/env python3
"""Summarize intervention MAD for metabolic vs autoimmune (RA/IBD) vs WBC."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results/tables"
METABOLIC = {
    "PGS000018",
    "PGS004696",
    "PGS004698",
    "PGS003897",
    "PGS002853",
    "PGS002858",
    "PGS003092",
    "PGS000014",
    "PGS004840",
}
AUTOIMMUNE = {"PGS004133", "PGS001288"}
WBC = {"PGS000191"}


def main() -> None:
    path = TABLES / "intervention_results.genomewide.csv"
    if not path.exists():
        raise SystemExit(f"missing {path}")
    df = pd.read_csv(path)
    mad = df[df["metric"] == "mean_abs_delta_EUR"].copy() if "metric" in df.columns else df.copy()
    if "pgs_id" not in mad.columns:
        raise SystemExit("pgs_id missing")
    val = next((c for c in ["value", "mad", "estimate"] if c in mad.columns), None)
    if val is None:
        raise SystemExit("value column missing")

    def clade(pid: str) -> str:
        if pid in WBC:
            return "wbc_duffy"
        if pid in AUTOIMMUNE:
            return "autoimmune_mhc"
        if pid in METABOLIC:
            return "metabolic"
        return "other"

    mad["clade"] = mad["pgs_id"].map(clade)
    g = mad.groupby(["clade", "mode"], as_index=False)[val].mean()
    out = TABLES / "intervention_mad_metabolic_vs_autoimmune.csv"
    g.to_csv(out, index=False, float_format="%.6g")
    print(g.to_string(index=False))
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
