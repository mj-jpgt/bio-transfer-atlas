"""
Rebuild slim PLINK2 .pvar files from local score_lookup.parquet caches.

Use when OneDrive leaves huge .pvar files as online-only stubs that PLINK cannot open.
Keeps the original stub as chrN.score.pvar.onedrive_stub if present.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data/interim/1000g_grch38"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Rebuild slim score.pvar from lookup parquet.")
    p.add_argument("--chroms", default="8-21")
    return p.parse_args()


def parse_chroms(spec: str) -> list[str]:
    if "-" in spec:
        a, b = spec.split("-", 1)
        return [str(x) for x in range(int(a), int(b) + 1)]
    return [str(int(x.strip())) for x in spec.split(",") if x.strip()]


def rebuild(chrom: str) -> None:
    lookup = INTERIM / f"chr{chrom}.score_lookup.parquet"
    pgen = INTERIM / f"chr{chrom}.score.pgen"
    pvar = INTERIM / f"chr{chrom}.score.pvar"
    if not lookup.exists():
        print(f"chr{chrom}: no lookup parquet — skip")
        return
    if not pgen.exists():
        print(f"chr{chrom}: no score.pgen — skip")
        return

    # If pvar already local and small enough (<500 MB), leave it
    if pvar.exists():
        try:
            attrs = subprocess.check_output(["attrib", str(pvar)], text=True)
            size = pvar.stat().st_size
            if "O" not in attrs and size < 500_000_000:
                print(f"chr{chrom}: local slim/ok pvar ({size/1e6:.0f} MB) — skip")
                return
        except Exception:
            pass
        stub = INTERIM / f"chr{chrom}.score.pvar.onedrive_stub"
        if not stub.exists() and pvar.stat().st_size > 500_000_000:
            pvar.rename(stub)
            print(f"chr{chrom}: moved large/stub pvar -> {stub.name}")
        elif pvar.exists():
            pvar.unlink()

    df = pd.read_parquet(lookup)
    print(f"chr{chrom}: writing slim pvar from {len(df):,} lookup rows ...", flush=True)
    with open(pvar, "w", encoding="utf-8", newline="\n") as f:
        f.write("#CHROM\tPOS\tID\tREF\tALT\n")
        for chr_, pos, ref, alt in zip(
            df["chr"].astype(str),
            df["pos"].astype(int),
            df["ref"].astype(str),
            df["alt"].astype(str),
        ):
            vid = f"{chr_}:{pos}:{ref}:{alt}"
            f.write(f"{chr_}\t{pos}\t{vid}\t{ref}\t{alt}\n")
    subprocess.run(["attrib", "-U", "+P", str(pvar)], check=False)
    print(f"chr{chrom}: wrote {pvar.name} ({pvar.stat().st_size/1e6:.1f} MB)")


def main() -> None:
    args = parse_args()
    for chrom in parse_chroms(args.chroms):
        rebuild(chrom)


if __name__ == "__main__":
    main()
