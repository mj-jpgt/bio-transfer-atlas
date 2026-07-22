"""
Download Ensembl GRCh38 release-110 chromosome GFF3 files for chr1-22.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data/raw/ensembl"
OUT_DIR.mkdir(parents=True, exist_ok=True)

URL_TEMPLATE = (
    "https://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens/"
    "Homo_sapiens.GRCh38.110.chromosome.{chrom}.gff3.gz"
)


def parse_chrom_list(spec: str) -> list[str]:
    s = spec.strip()
    if "-" in s:
        a, b = s.split("-", 1)
        return [str(x) for x in range(int(a), int(b) + 1)]
    return [str(int(x.strip())) for x in s.split(",") if x.strip()]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chroms", default="1-22", help="Chromosome range/list, e.g. 1-22 or 1,2,22")
    return p.parse_args()


def download_one(chrom: str) -> None:
    out = OUT_DIR / f"chr{chrom}.gff3.gz"
    if out.exists():
        print(f"chr{chrom}: cached -> {out.name}")
        return
    url = URL_TEMPLATE.format(chrom=chrom)
    from bta_curl import curl_bin

    cmd = [curl_bin(), "-fLsS", "-o", str(out), url]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Failed download for chr{chrom}: {r.stderr[-400:]}")
    print(f"chr{chrom}: downloaded -> {out.name}")


def main() -> None:
    args = parse_args()
    chroms = parse_chrom_list(args.chroms)
    print(f"Downloading Ensembl GFF3 chr files: {chroms}")
    for chrom in chroms:
        download_one(chrom)
    print("Done.")


if __name__ == "__main__":
    main()
