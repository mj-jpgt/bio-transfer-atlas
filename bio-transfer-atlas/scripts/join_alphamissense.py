"""
Stream AlphaMissense hg38 from GCS and join only to a small variant set (low RAM).

Never materializes the full ~643 MB TSV in pandas. Writes
data/annotations/alphamissense_grch38.parquet for matched variants.
"""
from __future__ import annotations

import argparse
import gzip
import io
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
ANN = ROOT / "data/annotations"
AM_URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--variant-list",
        default=str(ROOT / "data/modeling/_tmp_ldblock_associated_sample.parquet"),
        help="Parquet with variant_id; keep AM rows matching these only",
    )
    p.add_argument("--out", default=str(ANN / "alphamissense_grch38.parquet"))
    p.add_argument("--min-free-gb", type=float, default=1.5)
    return p.parse_args()


def free_gb() -> float:
    try:
        import psutil

        return psutil.virtual_memory().available / 1e9
    except Exception:
        pass
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / (1024**2)
    except OSError:
        pass
    return 99.0


def main() -> None:
    args = parse_args()
    avail = free_gb()
    print(f"RAM available: {avail:.2f} GB", flush=True)
    if avail < args.min_free_gb:
        raise SystemExit(
            f"Need >= {args.min_free_gb} GB free RAM (have {avail:.2f}). "
            "Close browsers/OneDrive sync and retry."
        )

    vpath = Path(args.variant_list)
    if not vpath.exists():
        raise SystemExit(f"Missing variant list {vpath}")
    wanted = set(pd.read_parquet(vpath, columns=["variant_id"])["variant_id"].astype(str).unique())
    print(f"Wanted variants: {len(wanted):,}", flush=True)

    print(f"Streaming {AM_URL} ...", flush=True)
    kept: list[dict] = []
    n_read = 0
    with requests.get(AM_URL, stream=True, timeout=120) as r:
        r.raise_for_status()
        # Decode gzip stream
        gz = gzip.GzipFile(fileobj=r.raw)
        text = io.TextIOWrapper(gz, encoding="utf-8", errors="replace")
        header = None
        for line in text:
            if not line.strip():
                continue
            if line.startswith("#"):
                if "CHROM" in line and "POS" in line:
                    header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if header is None:
                header = [
                    "CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "genome",
                    "uniprot_id",
                    "transcript_id",
                    "protein_variant",
                    "am_pathogenicity",
                    "am_class",
                ]
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            row = dict(zip(header, parts))
            chrom = str(row["CHROM"]).replace("chr", "")
            vid = f"{chrom}:{row['POS']}:{row['REF']}:{row['ALT']}"
            n_read += 1
            if vid in wanted:
                try:
                    score = float(row["am_pathogenicity"])
                except (TypeError, ValueError):
                    score = float("nan")
                kept.append(
                    {
                        "variant_id": vid,
                        "alphamissense_score": score,
                        "am_class": row.get("am_class", ""),
                    }
                )
            if n_read % 5_000_000 == 0:
                print(f"  scanned {n_read:,}; matched {len(kept):,}", flush=True)
                if free_gb() < 0.8:
                    print("  WARN low RAM — continuing carefully", flush=True)

    if not kept:
        raise SystemExit(f"No overlaps after scanning {n_read:,} AM rows")
    out = pd.DataFrame(kept).drop_duplicates("variant_id")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out, index=False)
    print(f"Saved {args.out} ({len(out):,} / {len(wanted):,} wanted; scanned {n_read:,})")
    print(out["alphamissense_score"].describe())


if __name__ == "__main__":
    main()
