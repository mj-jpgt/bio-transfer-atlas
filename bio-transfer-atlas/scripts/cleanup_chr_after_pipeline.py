"""
Delete regenerable per-chromosome caches after master table is written.

Keeps chr{N}.score.{pgen,pvar,psam}; removes LD .vcor, AF/LD/SEL feature files,
and interim QC/prune artifacts for that chromosome.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data/interim/1000g_grch38"
FEATURES = ROOT / "data/features"
CHR_RE = re.compile(r"chr(\d+)")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Free disk after a chromosome finishes pipeline.")
    p.add_argument("--chrom", required=True, help="Chromosome number, e.g. 8")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def chrom_in_name(name: str) -> int | None:
    m = CHR_RE.search(name)
    return int(m.group(1)) if m else None


def delete_matching(directory: Path, chrom: int, dry_run: bool) -> tuple[int, int]:
    if not directory.exists():
        return 0, 0
    n, b = 0, 0
    for f in list(directory.iterdir()):
        if not f.is_file():
            continue
        c = chrom_in_name(f.name)
        if c != chrom:
            continue
        b += f.stat().st_size
        if not dry_run:
            f.unlink()
        n += 1
    return n, b


def delete_ld_vcor(chrom: int, dry_run: bool) -> tuple[int, int]:
    ld_dir = FEATURES / "ld"
    if not ld_dir.exists():
        return 0, 0
    n, b = 0, 0
    for f in ld_dir.glob(f"*chr{chrom}*.vcor"):
        b += f.stat().st_size
        if not dry_run:
            f.unlink()
        n += 1
    return n, b


def delete_interim_extras(chrom: int, dry_run: bool) -> tuple[int, int]:
    if not INTERIM.exists():
        return 0, 0
    keep_prefixes = (
        f"chr{chrom}.score.pgen",
        f"chr{chrom}.score.pvar",
        f"chr{chrom}.score.psam",
    )
    n, b = 0, 0
    for f in list(INTERIM.iterdir()):
        if not f.is_file() or not f.name.startswith(f"chr{chrom}."):
            continue
        if f.name.startswith(keep_prefixes):
            continue
        b += f.stat().st_size
        if not dry_run:
            f.unlink()
        n += 1
    return n, b


def cleanup_chromosome(chrom: str, dry_run: bool = False) -> tuple[int, int]:
    c = int(chrom)
    total_n, total_b = 0, 0
    steps = [
        ("LD .vcor", *delete_ld_vcor(c, dry_run)),
        ("AF features", *delete_matching(FEATURES / "af", c, dry_run)),
        ("LD features", *delete_matching(FEATURES / "ld", c, dry_run)),
        ("SEL features", *delete_matching(FEATURES / "selection", c, dry_run)),
        ("Interim extras (keep score.pgen)", *delete_interim_extras(c, dry_run)),
    ]
    label = "DRY RUN" if dry_run else "cleanup"
    print(f"chr{c} {label}:")
    for name, n, b in steps:
        total_n += n
        total_b += b
        print(f"  {name}: {n:,} files, {b / (1024**3):.2f} GB")
    print(f"  total: {total_n:,} files, {total_b / (1024**3):.2f} GB")
    return total_n, total_b


def main() -> None:
    args = parse_args()
    chrom = str(int(args.chrom))
    master = ROOT / f"data/modeling/master_variant_table.chr{chrom}.parquet"
    if not master.exists() and not args.dry_run:
        raise SystemExit(f"Refusing cleanup: missing {master}")
    cleanup_chromosome(chrom, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
