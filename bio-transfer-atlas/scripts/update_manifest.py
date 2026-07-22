"""
Scan data/raw/ and update data/raw/MANIFEST.tsv with provenance info.

Usage:
    python scripts/update_manifest.py
"""

import csv
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

from loguru import logger

MANIFEST_COLS = [
    "source",
    "dataset",
    "file",
    "url",
    "date_downloaded",
    "sha256",
    "license",
    "notes",
]

SOURCE_MAP = {
    "1000g": ("IGSR / 1000 Genomes", "public domain"),
    "pgs_catalog": ("PGS Catalog", "see score-level license in catalog"),
    "reactome": ("Reactome", "Creative Commons Attribution 4.0"),
    "panukbb": ("Pan-UKBB", "CC BY 4.0"),
    "bbj": ("BioBank Japan PheWeb", "see BBJ terms"),
    "finngen": ("FinnGen", "see FinnGen terms"),
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    raw_dir = root / "data" / "raw"
    manifest_path = raw_dir / "MANIFEST.tsv"

    existing: dict[str, dict] = {}
    if manifest_path.exists():
        with open(manifest_path, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                existing[row["file"]] = row

    rows = []
    for fpath in sorted(raw_dir.rglob("*")):
        if not fpath.is_file():
            continue
        rel = fpath.relative_to(raw_dir).as_posix()
        if rel == "MANIFEST.tsv":
            continue

        source_key = rel.split("/")[0]
        source, license_ = SOURCE_MAP.get(source_key, ("unknown", "unknown"))

        if rel in existing:
            row = existing[rel]
            if not row.get("sha256"):
                row["sha256"] = sha256_file(fpath)
        else:
            stat = fpath.stat()
            row = {
                "source": source,
                "dataset": source_key,
                "file": rel,
                "url": "",
                "date_downloaded": datetime.fromtimestamp(
                    stat.st_mtime, tz=timezone.utc
                ).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "sha256": sha256_file(fpath),
                "license": license_,
                "notes": "",
            }
        rows.append(row)

    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_COLS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    logger.success(f"MANIFEST updated: {manifest_path}  ({len(rows)} files)")


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="INFO")
    main()
