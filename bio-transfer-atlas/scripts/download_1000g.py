"""
Download 1000 Genomes Phase 3 data needed for bio-transfer-atlas.

Downloads:
  - Sample metadata panel file
  - chr22 VCF (start small; set CHROMOSOMES env var to extend)
  - VCF index (.tbi)

Usage:
    python scripts/download_1000g.py
    CHROMOSOMES=1,2,22 python scripts/download_1000g.py
"""

import hashlib
import os
import sys
import time
from pathlib import Path

import requests
from loguru import logger
from tqdm import tqdm

BASE_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

METADATA_FILE = "integrated_call_samples_v3.20130502.ALL.panel"

VCF_TEMPLATE = (
    "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)

CHROMOSOMES = [
    int(c.strip())
    for c in os.environ.get("CHROMOSOMES", "22").split(",")
    if c.strip()
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def download_file(url: str, dest: Path, desc: str = "") -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        logger.info(f"Already exists, skipping: {dest}")
        return

    logger.info(f"Downloading {desc or url} -> {dest}")
    r = requests.get(url, stream=True, timeout=60)
    r.raise_for_status()

    total = int(r.headers.get("content-length", 0))
    with open(dest, "wb") as f, tqdm(
        total=total, unit="B", unit_scale=True, desc=dest.name
    ) as bar:
        for chunk in r.iter_content(chunk_size=131072):
            f.write(chunk)
            bar.update(len(chunk))

    logger.success(f"Saved: {dest}  sha256={sha256_file(dest)[:12]}...")


def main() -> None:
    root = Path(__file__).resolve().parents[1]

    metadata_dir = root / "data" / "raw" / "1000g" / "metadata"
    vcf_dir = root / "data" / "raw" / "1000g" / "vcf"

    metadata_url = f"{BASE_URL}/{METADATA_FILE}"
    download_file(metadata_url, metadata_dir / METADATA_FILE, "1000G sample metadata")

    for chrom in CHROMOSOMES:
        vcf_name = VCF_TEMPLATE.format(chrom=chrom)
        vcf_url = f"{BASE_URL}/{vcf_name}"
        tbi_url = f"{vcf_url}.tbi"

        download_file(vcf_url, vcf_dir / vcf_name, f"1000G chr{chrom} VCF")
        download_file(tbi_url, vcf_dir / f"{vcf_name}.tbi", f"1000G chr{chrom} TBI")

    logger.success(f"1000G downloads complete for chromosomes: {CHROMOSOMES}")


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="INFO")
    main()
