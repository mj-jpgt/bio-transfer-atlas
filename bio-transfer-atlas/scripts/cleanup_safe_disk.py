"""
Delete pipeline artifacts for chromosomes with finished master tables (1-7, 22).
Does NOT touch chr8-21 (pipeline still pending).
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DONE = {1, 2, 3, 4, 5, 6, 7, 22}
INTERIM_DELETE = {1, 2, 3, 4, 5, 6, 7}  # chr22: delete extras only, keep score.pgen
INTERIM = ROOT / "data/interim/1000g_grch38"
FEATURES = ROOT / "data/features"

CHR_RE = re.compile(r"chr(\d+)")


def chrom_in_name(name: str) -> int | None:
    m = CHR_RE.search(name)
    return int(m.group(1)) if m else None


def is_done_chr_name(name: str) -> bool:
    c = chrom_in_name(name)
    return c is not None and c in DONE


def delete_glob(directory: Path, pattern: str) -> tuple[int, int]:
    if not directory.exists():
        return 0, 0
    n, b = 0, 0
    for f in directory.glob(pattern):
        if f.is_file():
            b += f.stat().st_size
            f.unlink()
            n += 1
    return n, b


def delete_features_done_chrs() -> tuple[int, int]:
    n, b = 0, 0
    for sub in ("af", "ld", "selection"):
        d = FEATURES / sub
        if not d.exists():
            continue
        for f in list(d.iterdir()):
            if f.is_file() and is_done_chr_name(f.name):
                b += f.stat().st_size
                f.unlink()
                n += 1
    return n, b


def delete_interim_done_chrs() -> tuple[int, int]:
    n, b = 0, 0
    if not INTERIM.exists():
        return 0, 0
    for f in list(INTERIM.iterdir()):
        if not f.is_file():
            continue
        name = f.name
        c = chrom_in_name(name)
        if c is not None and c in INTERIM_DELETE and name.startswith("chr"):
            b += f.stat().st_size
            f.unlink()
            n += 1
    return n, b


def delete_chr22_extras() -> tuple[int, int]:
    n, b = 0, 0
    keep_prefixes = ("chr22.score.pgen", "chr22.score.pvar", "chr22.score.psam")
    for f in list(INTERIM.iterdir()):
        if not f.is_file() or not f.name.startswith("chr22."):
            continue
        if f.name.startswith(keep_prefixes):
            continue
        b += f.stat().st_size
        f.unlink()
        n += 1
    return n, b


def delete_interim_junk() -> tuple[int, int]:
    n, b = 0, 0
    if not INTERIM.exists():
        return 0, 0
    patterns = [
        "genome_wide_grch38.*",
        "all_chrs_grch38.*",
        "pmerge_list_grch38.txt",
        "_*",
    ]
    for pat in patterns:
        for f in INTERIM.glob(pat):
            if f.is_file():
                b += f.stat().st_size
                f.unlink()
                n += 1
    return n, b


def main() -> None:
    total_n, total_b = 0, 0
    steps = [
        ("LD .vcor caches (all)", *delete_glob(FEATURES / "ld", "*.vcor")),
        ("LD .log files (all)", *delete_glob(FEATURES / "ld", "*.log")),
        ("AF .log files (all)", *delete_glob(FEATURES / "af", "*.log")),
        ("Features chr1-7,22", *delete_features_done_chrs()),
        ("Interim chr1-7", *delete_interim_done_chrs()),
        ("Interim chr22 extras (keep score.pgen)", *delete_chr22_extras()),
        ("Interim genome-wide / test junk", *delete_interim_junk()),
    ]
    print("Safe cleanup (masters done: chr 1-7, 22; chr8-21 untouched)\n")
    for label, n, b in steps:
        total_n += n
        total_b += b
        print(f"  {label}: {n:,} files, {b / (1024**3):.2f} GB")
    print(f"\nTotal removed: {total_n:,} files, {total_b / (1024**3):.2f} GB")


if __name__ == "__main__":
    main()
