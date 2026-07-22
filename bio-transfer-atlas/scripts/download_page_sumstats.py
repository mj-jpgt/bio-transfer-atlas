"""
Download PAGE (Wojcik 2019) GWAS Catalog sumstats, keep one chromosome, write parquet
compatible with run_external_sumstat_validation.py naming: page_{trait}_{anc}.parquet
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/raw/external_sumstats"

# Harmonised build37 files (smaller than .h); streamed and chrom-filtered.
PAGE_FILES = {
    "LDL": {
        "url": (
            "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            "GCST008001-GCST009000/GCST008037/harmonised/"
            "31217584-GCST008037-EFO_0004611-build37.f.tsv.gz"
        ),
        "anc_label": "AFR",  # PAGE multi-ancestry non-EUR cohort; label for pairing
        "note": "PAGE Wojcik2019 LDL (multi-ancestry; labeled AFR for EUR-vs-PAGE pair)",
    },
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--trait", default="LDL", choices=list(PAGE_FILES))
    p.add_argument("--chrom", default="22")
    p.add_argument("--max-rows", type=int, default=0, help="Optional cap after chrom filter")
    return p.parse_args()


def normalize_chr(v) -> str:
    s = str(v).strip()
    return s[3:] if s.lower().startswith("chr") else s


def main() -> None:
    args = parse_args()
    meta = PAGE_FILES[args.trait]
    OUT.mkdir(parents=True, exist_ok=True)
    out_pq = OUT / f"page_{args.trait}_{meta['anc_label']}.parquet"
    if out_pq.exists() and out_pq.stat().st_size > 100_000:
        print(f"cached {out_pq}")
        return

    print(f"Streaming {meta['url']}", flush=True)
    with requests.get(meta["url"], stream=True, timeout=120) as r:
        r.raise_for_status()
        # Write to temp then stream-read with pandas
        tmp = OUT / f"_tmp_page_{args.trait}.tsv.gz"
        n = 0
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                if chunk:
                    f.write(chunk)
                    n += len(chunk)
                    if n % (64 * 1024 * 1024) < 8 * 1024 * 1024:
                        print(f"  downloaded {n / 1e6:.0f} MB ...", flush=True)
        print(f"  download complete {n / 1e6:.1f} MB", flush=True)

    # Peek header
    with gzip.open(tmp, "rt", encoding="utf-8", errors="replace") as gz:
        header = gz.readline().rstrip("\n").split("\t")
    print(f"  columns ({len(header)}): {header[:12]} ...", flush=True)

    # Map common GWAS Catalog harmonised names
    lower = {c.lower(): c for c in header}

    def pick(*cands: str) -> str | None:
        for c in cands:
            if c.lower() in lower:
                return lower[c.lower()]
        return None

    chrom_c = pick("chromosome", "chr", "hm_chrom", "chr_id")
    pos_c = pick("base_pair_location", "pos", "hm_pos", "bp", "position")
    ea_c = pick("effect_allele", "hm_effect_allele", "a1", "alt")
    oa_c = pick("other_allele", "hm_other_allele", "a2", "ref")
    beta_c = pick("beta", "hm_beta", "effect", "effect_weight")
    se_c = pick("standard_error", "hm_se", "se")
    p_c = pick("p_value", "hm_pvalue", "pval", "p")

    need = [c for c in [chrom_c, pos_c, ea_c, oa_c, beta_c, se_c, p_c] if c]
    chunks = []
    n_keep = 0
    for chunk in pd.read_csv(
        tmp,
        sep="\t",
        compression="gzip",
        usecols=need,
        chunksize=200_000,
        dtype=str,
        low_memory=False,
    ):
        chunk = chunk[chunk[chrom_c].map(normalize_chr) == str(args.chrom)]
        if chunk.empty:
            continue
        chunks.append(chunk)
        n_keep += len(chunk)
        if args.max_rows and n_keep >= args.max_rows:
            break
        if n_keep and n_keep % 500_000 < 200_000:
            print(f"  kept {n_keep:,} chr{args.chrom} rows", flush=True)

    if not chunks:
        raise SystemExit(f"No chr{args.chrom} rows found")
    df = pd.concat(chunks, ignore_index=True)
    if args.max_rows:
        df = df.head(args.max_rows)

    # Build GRCh38-style variant_id if possible; PAGE build37 — use chr:pos:ref:alt as-is
    # and also emit a lift-free join key on chrom:pos for Pan-UKB matching via pos37 if needed.
    df["chrom"] = df[chrom_c].map(normalize_chr)
    df["pos"] = pd.to_numeric(df[pos_c], errors="coerce").astype("Int64")
    df["effect_allele"] = df[ea_c].astype(str).str.upper()
    df["other_allele"] = df[oa_c].astype(str).str.upper()
    df["beta"] = pd.to_numeric(df[beta_c], errors="coerce")
    df["se"] = pd.to_numeric(df[se_c], errors="coerce") if se_c else float("nan")
    df["pval"] = pd.to_numeric(df[p_c], errors="coerce") if p_c else float("nan")
    df = df.dropna(subset=["pos", "beta"])
    df["variant_id"] = (
        df["chrom"].astype(str)
        + ":"
        + df["pos"].astype(str)
        + ":"
        + df["other_allele"]
        + ":"
        + df["effect_allele"]
    )
    out = df[["variant_id", "beta", "se", "pval", "chrom", "pos"]].copy()
    out.to_parquet(out_pq, index=False)
    note = OUT / f"page_{args.trait}_{meta['anc_label']}.note.txt"
    note.write_text(meta["note"] + f"\nchrom={args.chrom}\nn={len(out)}\n", encoding="utf-8")
    tmp.unlink(missing_ok=True)
    print(f"Saved {out_pq} ({len(out):,} rows)")


if __name__ == "__main__":
    main()
