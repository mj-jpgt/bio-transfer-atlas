"""
Compute gene-constraint features per variant using gnomAD v4.1.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/reference/gnomad.v4.1.constraint_metrics.tsv"
URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"

OUT_DIR = ROOT / "data/features/selection"
OUT_DIR.mkdir(parents=True, exist_ok=True)
RAW.parent.mkdir(parents=True, exist_ok=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22", help="Chromosome number/name (e.g. 22 or chr22)")
    return p.parse_args()


def normalize_chr(value: str) -> str:
    s = str(value).strip()
    return s[3:] if s.lower().startswith("chr") else s


def resolve_v2g_path(chrom: str) -> Path:
    p = ROOT / "data/annotations" / f"variant_to_gene.chr{chrom}.parquet"
    if p.exists():
        return p
    if chrom == "22":
        legacy = ROOT / "data/annotations/variant_to_gene.parquet"
        if legacy.exists():
            return legacy
    raise FileNotFoundError(f"Missing variant_to_gene for chr{chrom}: {p}")


def maybe_write_legacy(chrom: str, src: Path) -> None:
    if chrom == "22":
        pd.read_parquet(src).to_parquet(OUT_DIR / "constraint_features.parquet", index=False)


def as_bool(x):
    if pd.isna(x):
        return False
    s = str(x).strip().lower()
    return s in {"true", "1", "t", "yes"}


if __name__ == "__main__":
    args = parse_args()
    chrom = normalize_chr(args.chrom)
    v2g_path = resolve_v2g_path(chrom)
    out_path = OUT_DIR / f"constraint_features.chr{chrom}.parquet"

    if not RAW.exists():
        print(f"Downloading gnomAD constraint table -> {RAW}")
        pd.read_csv(URL, sep="\t", dtype=str).to_csv(RAW, sep="\t", index=False)

    usecols = [
        "gene_id",
        "transcript",
        "canonical",
        "mane_select",
        "lof.oe_ci.upper",
        "lof.pLI",
        "mis.z_score",
    ]
    c = pd.read_csv(RAW, sep="\t", dtype=str, usecols=usecols)
    c = c.rename(
        columns={
            "gene_id": "ensg",
            "lof.oe_ci.upper": "LOEUF",
            "lof.pLI": "pLI",
            "mis.z_score": "mis_z",
        }
    )
    c["LOEUF"] = pd.to_numeric(c["LOEUF"], errors="coerce")
    c["pLI"] = pd.to_numeric(c["pLI"], errors="coerce")
    c["mis_z"] = pd.to_numeric(c["mis_z"], errors="coerce")
    c["mane_select"] = c["mane_select"].map(as_bool)
    c["canonical"] = c["canonical"].map(as_bool)

    # transcript priority: MANE > canonical > lowest LOEUF as fallback
    c["priority"] = 2
    c.loc[c["canonical"], "priority"] = 1
    c.loc[c["mane_select"], "priority"] = 0
    c = c.sort_values(["ensg", "priority", "LOEUF"], na_position="last")
    gene = c.groupby("ensg", as_index=False).first()[["ensg", "LOEUF", "pLI", "mis_z", "transcript"]]

    v2g = pd.read_parquet(v2g_path)[["variant_id", "ensg"]].drop_duplicates()
    m = v2g.merge(gene, on="ensg", how="left")

    # For variants mapping to multiple genes, pick most constrained (min LOEUF).
    m = m.sort_values(["variant_id", "LOEUF"], na_position="last")
    best = m.groupby("variant_id", as_index=False).first()
    out = best.rename(columns={"transcript": "constraint_transcript"})
    out = out[["variant_id", "LOEUF", "pLI", "mis_z", "constraint_transcript"]]
    out.to_parquet(out_path, index=False)
    maybe_write_legacy(chrom, out_path)

    print(f"Saved {out_path} ({len(out):,} variants)")
    print(f"Coverage LOEUF non-null: {out['LOEUF'].notna().mean()*100:.1f}%")
