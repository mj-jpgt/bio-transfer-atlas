"""
Download BioBank Japan GWAS and keep chr22 only.

This script handles mixed BBJ formats:
- v3 zip bundles containing an inner .txt.gz member
- v5 direct .txt.gz files
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Dict, Iterable, Iterator, List

import pandas as pd

from gwas_utils import finalize_gwas_frame, normalize_chr, read_tsv_gz_chunks, read_zip_member_chunks

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22", help="Chromosome number/name (e.g. 22 or chr22)")
    return p.parse_args()

# Stable links from https://pheweb.jp/downloads rows (BBJ cohort).
TRAIT_URLS: Dict[str, str] = {
    "T2D": "https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.T2D.v1.zip",
    "CAD": "https://humandbs.dbcls.jp/files/hum0197/hum0197.v5.gwas.CAD.v1.txt.gz",
    "BMI": "https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.BMI.v1.zip",
    "LDL": "https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.LDLC.v1.zip",
}

ALIASES: Dict[str, List[str]] = {
    "chr": ["CHR", "#chrom", "#CHROM", "chrom", "chr"],
    "pos": ["POS", "pos", "BP", "bp"],
    "ref": ["Allele1", "ALLELE1", "REF", "ref", "A1", "allele1"],
    "alt": ["Allele2", "ALLELE2", "ALLELE0", "ALT", "alt", "A2", "allele2"],
    "beta": ["BETA", "beta", "Effect", "effect", "B"],
    "se": ["SE", "se", "SEBETA", "sebeta"],
    "pval": ["p.value", "P", "pval", "PVALUE", "p.value.NA", "P_BOLT_LMM_INF", "P_LINREG"],
    "af": ["AF_Allele2", "A1FREQ", "af_alt", "AF_ALT", "Frq", "EAF", "AF"],
}


def resolve_aliases(columns: Iterable[str]) -> Dict[str, str]:
    cols = list(columns)
    out: Dict[str, str] = {}
    for canonical, choices in ALIASES.items():
        hit = next((c for c in choices if c in cols), None)
        if hit is not None:
            out[canonical] = hit
    required = ["chr", "pos", "ref", "alt", "beta", "se", "pval"]
    missing = [c for c in required if c not in out]
    if missing:
        raise RuntimeError(f"Could not resolve required columns {missing}; saw columns: {cols[:30]}")
    if "af" not in out:
        out["af"] = out["pval"]  # temporary placeholder; overwritten to NaN later
    return out


def iter_chunks(url: str) -> Iterator[pd.DataFrame]:
    if url.endswith(".zip"):
        yield from read_zip_member_chunks(url, chunksize=250_000)
    else:
        yield from read_tsv_gz_chunks(url, chunksize=350_000)


def process_trait(trait: str, url: str, chrom: str = "22") -> Path:
    chrom = normalize_chr(chrom)
    out_dir = ROOT / "data/raw/bbj" / f"chr{chrom}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{trait}.chr{chrom}.parquet"
    if out_path.exists():
        print(f"{trait}: cached -> {out_path.name}")
        return out_path

    print(f"{trait}: reading {url}")
    frames: List[pd.DataFrame] = []
    alias_map: Dict[str, str] | None = None
    max_attempts = 4
    for attempt in range(1, max_attempts + 1):
        try:
            frames = []
            alias_map = None
            for i, chunk in enumerate(iter_chunks(url), start=1):
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
            break
        except Exception as e:
            msg = str(e)
            transient = (
                isinstance(e, (ConnectionAbortedError, ConnectionResetError, TimeoutError, OSError))
                or "10053" in msg
                or "Connection aborted" in msg
                or "timed out" in msg.lower()
            )
            if (not transient) or attempt == max_attempts:
                raise
            wait_s = 5 * attempt
            print(
                f"{trait}: transient network/read failure on attempt {attempt}/{max_attempts}, "
                f"retrying in {wait_s}s ... ({e})"
            )
            time.sleep(wait_s)

    if not frames:
        raise RuntimeError(f"{trait}: no chr{chrom} rows found")
    df = pd.concat(frames, ignore_index=True)
    df = finalize_gwas_frame(df)
    if "af" in df and alias_map and alias_map.get("af") == alias_map.get("pval"):
        df["af"] = pd.NA
    df = df[["chr", "pos", "ref", "alt", "beta", "se", "pval", "af"]]
    df.to_parquet(out_path, index=False)
    print(f"{trait}: saved {len(df):,} rows -> {out_path.name}")
    return out_path


if __name__ == "__main__":
    args = parse_args()
    for trait, url in TRAIT_URLS.items():
        process_trait(trait, url, chrom=args.chrom)
    print("Done.")
