"""
Pre-flight checks before long genome-wide runs (RAM, score pgens, OneDrive note).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from rebuild_score_pfile_conservative import free_ram_gb, score_pfile_ok
from intervention_common import available_score_chroms, GW_CHROMS

INTERIM = ROOT / "data/interim/1000g_grch38"
VCF_DIR = ROOT / "data/raw/1000g/vcf_grch38"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Genome-wide preflight checks.")
    p.add_argument("--min-free-gb", type=float, default=4.0)
    p.add_argument("--rebuild-chr22", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    free = free_ram_gb()
    print(f"Free RAM: {free:.1f} GB (recommend >= {args.min_free_gb} GB)")
    if free < args.min_free_gb:
        print("WARNING: Low RAM — close other apps before PLINK/LD steps.")

    score_chroms = available_score_chroms()
    missing = [c for c in GW_CHROMS if c not in score_chroms]
    print(f"Score pgens present: {len(score_chroms)}/22")
    if missing:
        print(f"  Missing score pgens: {', '.join(missing)}")

    if args.rebuild_chr22 and "22" in missing:
        vcf = VCF_DIR / "1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        if not vcf.exists():
            print(f"chr22 VCF missing locally: {vcf}")
            print("  Restore from Hugging Face (mj0jpgg/fairgen) or IGSR before rebuild.")
        else:
            import subprocess
            py = sys.executable
            cmd = [
                py, "scripts/rebuild_score_pfile_conservative.py",
                "--chrom", "22", "--memory-mb", "2048", "--threads", "1",
                "--min-free-gb", str(args.min_free_gb),
            ]
            print("Rebuilding chr22.score ...")
            subprocess.run(cmd, cwd=ROOT, check=False)

    if score_pfile_ok("22", plink_check=False):
        print("chr22.score: OK (file check)")
    elif (INTERIM / "chr22.score.pgen").exists():
        print("chr22.score: present (PLINK sanity skipped)")

    print("\nOneDrive: pause sync on data/interim/ during PLINK to avoid cloud-stub errors.")
    print("Run pipeline: python scripts/run_genomewide_conservative.py --chroms 8-21 --skip-score-rebuild ...")


if __name__ == "__main__":
    main()
