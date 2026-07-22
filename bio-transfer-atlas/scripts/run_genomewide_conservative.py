"""
Low-RAM, resumable genome-wide pipeline driver.

- Skips chromosomes with an existing master table
- Rebuilds broken score pgens on local disk before AF/LD
- Requires minimum free RAM before each chromosome
- No GWAS re-downloads; single-threaded PLINK where possible
"""
from __future__ import annotations

import argparse
import gc
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def parse_chrom_list(spec: str) -> list[str]:
    s = spec.strip()
    if "-" in s:
        a, b = s.split("-", 1)
        return [str(x) for x in range(int(a), int(b) + 1)]
    return [str(int(x.strip())) for x in s.split(",") if x.strip()]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Conservative genome-wide pipeline.")
    p.add_argument("--chroms", default="7-21")
    p.add_argument("--min-free-gb", type=float, default=4.0)
    p.add_argument("--plink-memory-mb", type=int, default=2048)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--skip-pathways", action="store_true")
    p.add_argument("--skip-harmonize", action="store_true")
    p.add_argument("--allow-qc-fallback", action="store_true",
                   help="Use qc.id clone if VCF import fails for score pfile")
    p.add_argument("--qc-fallback-only", action="store_true",
                   help="Score pfile repair: qc.id clone only (low RAM)")
    p.add_argument("--skip-score-rebuild", action="store_true",
                   help="Skip score pfile repair if .pgen/.pvar/.psam exist")
    p.add_argument("--python", default=sys.executable)
    return p.parse_args()


def score_pfile_exists(chrom: str) -> bool:
    interim = ROOT / "data/interim/1000g_grch38"
    for ext in (".pgen", ".pvar", ".psam"):
        p = interim / f"chr{chrom}.score{ext}"
        if not p.exists() or p.stat().st_size < 1000:
            return False
    return True


def run(cmd: list[str], desc: str) -> None:
    print(f"\n[{desc}]", flush=True)
    print("  " + " ".join(cmd), flush=True)
    env = {**__import__("os").environ, "PYTHONUNBUFFERED": "1"}
    r = subprocess.run(cmd, cwd=ROOT, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"FAILED ({r.returncode}): {desc}")


def master_done(chrom: str) -> bool:
    return (ROOT / f"data/modeling/master_variant_table.chr{chrom}.parquet").exists()


def ensure_score_pfile(chrom: str, py: str, memory_mb: int, threads: int, min_gb: float,
                       qc_fb: bool, qc_only: bool) -> None:
    cmd = [
        py, "scripts/rebuild_score_pfile_conservative.py",
        "--chrom", chrom,
        "--memory-mb", str(memory_mb),
        "--threads", str(threads),
        "--min-free-gb", str(min_gb),
    ]
    if qc_only:
        cmd.append("--qc-fallback-only")
    elif qc_fb:
        cmd.append("--allow-qc-fallback")
    run(cmd, f"Ensure score pfile chr{chrom}")


def main() -> None:
    args = parse_args()
    chroms = parse_chrom_list(args.chroms)
    py = args.python
    log_dir = ROOT / "logs"
    log_dir.mkdir(exist_ok=True)

    print(f"Conservative genome-wide run: {chroms}", flush=True)
    print(f"  min_free_gb={args.min_free_gb}  plink_memory={args.plink_memory_mb}  threads={args.threads}", flush=True)

    if not args.skip_harmonize:
        run([py, "scripts/harmonize_pgs_genomewide.py"], "Harmonize PGS (genome-wide)")

    for chrom in chroms:
        if master_done(chrom):
            print(f"\nchr{chrom}: master table exists — skip", flush=True)
            continue

        print(f"\n{'='*60}\nchr{chrom} START\n{'='*60}", flush=True)
        if args.skip_score_rebuild and score_pfile_exists(chrom):
            print(f"chr{chrom}: score pfile present — skip rebuild", flush=True)
        else:
            ensure_score_pfile(
                chrom, py, args.plink_memory_mb, args.threads,
                args.min_free_gb, args.allow_qc_fallback, args.qc_fallback_only,
            )

        run([py, "scripts/compute_af_features.py", "--chrom", chrom, "--memory-mb", str(args.plink_memory_mb)],
            f"AF chr{chrom}")
        gc.collect()
        run([py, "scripts/compute_ld_features.py", "--chrom", chrom, "--threads", str(args.threads),
             "--memory-mb", str(args.plink_memory_mb)], f"LD chr{chrom}")
        gc.collect()
        run([py, "scripts/compute_selection_features.py", "--chrom", chrom], f"SEL chr{chrom}")
        run([py, "scripts/build_concordance_labels_multisource.py", "--chrom", chrom], f"Labels chr{chrom}")
        run([py, "scripts/build_master_table.py", "--chrom", chrom], f"Master chr{chrom}")
        if not args.skip_pathways:
            run([py, "scripts/build_pathway_risk.py", "--chrom", chrom], f"Pathway chr{chrom}")
            run([py, "scripts/compute_constraint_features.py", "--chrom", chrom], f"Constraint chr{chrom}")
            run([py, "scripts/build_master_table.py", "--chrom", chrom], f"Master refresh chr{chrom}")
        gc.collect()
        if master_done(chrom):
            run([py, "scripts/cleanup_chr_after_pipeline.py", "--chrom", chrom],
                f"Disk cleanup chr{chrom}")
        print(f"chr{chrom} DONE", flush=True)

    print("\nPer-chromosome steps complete. Run concat/eval separately when all chrs ready:", flush=True)
    print("  python scripts/concat_genomewide.py ...", flush=True)
    print("  python scripts/run_partial_eval.py --chroms ...", flush=True)


if __name__ == "__main__":
    main()
