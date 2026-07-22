"""
Download Reactome gene-pathway mapping files for bio-transfer-atlas.

Downloads:
  - NCBI2Reactome.txt         -- Entrez gene ID -> Reactome pathway
  - Ensembl2Reactome.txt      -- Ensembl gene ID -> Reactome pathway (curated)
  - Ensembl2Reactome_All_Levels.txt  -- includes all hierarchy levels
  - UniProt2Reactome.txt      -- UniProt -> Reactome pathway

Usage:
    python scripts/download_reactome.py
"""

import hashlib
import sys
from pathlib import Path

import requests
from loguru import logger
from tqdm import tqdm

REACTOME_FILES = [
    "NCBI2Reactome.txt",
    "Ensembl2Reactome.txt",
    "Ensembl2Reactome_All_Levels.txt",
    "UniProt2Reactome.txt",
]

BASE_URL = "https://reactome.org/download/current"
VERSION_URL = "https://reactome.org/ContentService/data/database/version"


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
    r = requests.get(url, stream=True, timeout=120)
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
    dest_dir = root / "data" / "raw" / "reactome"
    dest_dir.mkdir(parents=True, exist_ok=True)

    version_file = dest_dir / "version.txt"
    if not version_file.exists():
        r = requests.get(VERSION_URL, timeout=30)
        r.raise_for_status()
        version_file.write_text(r.text.strip())
        logger.info(f"Reactome version: {r.text.strip()}")
    else:
        logger.info(f"Reactome version cached: {version_file.read_text().strip()}")

    for fname in REACTOME_FILES:
        url = f"{BASE_URL}/{fname}"
        download_file(url, dest_dir / fname, fname)

    logger.success("Reactome downloads complete.")


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="INFO")
    main()
