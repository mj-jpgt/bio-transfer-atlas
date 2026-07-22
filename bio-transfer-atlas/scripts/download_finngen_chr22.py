"""
Download FinnGen summary statistics and keep chr22 only.

Notes:
- FinnGen R12 manifest is public and used by default.
- R13 may require form-based access instructions depending on endpoint/bucket.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd

from gwas_utils import finalize_gwas_frame, normalize_chr, read_tsv_gz_chunks

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22", help="Chromosome number/name (e.g. 22 or chr22)")
    return p.parse_args()

MANIFEST_URL = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv"

# FinnGen disease-focused endpoints; BMI/LDL are handled as meta sources.
TRAIT_ENDPOINT = {
    "T2D": "T2D",
    "CAD": "I9_CHD",
}

ALIASES: Dict[str, List[str]] = {
    "chr": ["#chrom", "chr", "CHROM", "#CHR"],
    "pos": ["pos", "POS"],
    "ref": ["ref", "REF", "allele1"],
    "alt": ["alt", "ALT", "allele2"],
    "beta": ["beta", "BETA"],
    "sebeta": ["sebeta", "SEBETA", "se"],
    "pval": ["pval", "P", "p.value"],
    "af_alt": ["af_alt", "AF_ALT", "af"],
}


def resolve_aliases(columns: Iterable[str]) -> Dict[str, str]:
    cols = list(columns)
    out: Dict[str, str] = {}
    for canonical, choices in ALIASES.items():
        hit = next((c for c in choices if c in cols), None)
        if hit is not None:
            out[canonical] = hit
    missing = [c for c in ["chr", "pos", "ref", "alt", "beta", "sebeta", "pval"] if c not in out]
    if missing:
        raise RuntimeError(f"Missing FinnGen columns {missing}; saw {cols[:30]}")
    if "af_alt" not in out:
        out["af_alt"] = out["pval"]  # temporary placeholder
    return out


def process_trait(trait: str, url: str, chrom: str = "22") -> Path:
    chrom = normalize_chr(chrom)
    out_dir = ROOT / "data/raw/finngen" / f"chr{chrom}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{trait}.chr{chrom}.parquet"
    if out_path.exists():
        print(f"{trait}: cached -> {out_path.name}")
        return out_path

    frames = []
    alias_map = None
    for i, chunk in enumerate(read_tsv_gz_chunks(url, chunksize=450_000), start=1):
        if alias_map is None:
            alias_map = resolve_aliases(chunk.columns)
            print(f"{trait}: resolved columns -> {alias_map}")
        keep = [alias_map[k] for k in alias_map]
        c = chunk[keep].rename(columns={v: k for k, v in alias_map.items()})
        c = c[c["chr"].map(normalize_chr) == str(chrom)]
        if len(c):
            frames.append(c)
        if i % 10 == 0:
            print(f"{trait}: processed {i} chunks")
    if not frames:
        raise RuntimeError(f"{trait}: no chr{chrom} rows found")
    df = pd.concat(frames, ignore_index=True)
    df = finalize_gwas_frame(df)
    if alias_map and alias_map.get("af_alt") == alias_map.get("pval"):
        df["af_alt"] = pd.NA
    df["se"] = df["sebeta"]
    df = df[["chr", "pos", "ref", "alt", "beta", "se", "sebeta", "pval", "af_alt"]]
    df.to_parquet(out_path, index=False)
    print(f"{trait}: saved {len(df):,} rows -> {out_path.name}")
    return out_path


if __name__ == "__main__":
    args = parse_args()
    man = pd.read_csv(MANIFEST_URL, sep="\t", dtype=str)
    for trait, endpoint in TRAIT_ENDPOINT.items():
        sub = man[man["phenocode"] == endpoint]
        if sub.empty:
            print(f"{trait}: endpoint {endpoint} not found in manifest; skipping")
            continue
        url = sub.iloc[0]["path_https"]
        process_trait(trait, url, chrom=args.chrom)
    print("Done.")
