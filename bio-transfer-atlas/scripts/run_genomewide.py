"""
Run the chr-by-chr portability pipeline and concatenate genome-wide artifacts.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SEED = 719


def parse_chrom_list(spec: str) -> list[str]:
    s = spec.strip()
    if "-" in s:
        a, b = s.split("-", 1)
        return [str(x) for x in range(int(a), int(b) + 1)]
    return [str(int(x.strip())) for x in s.split(",") if x.strip()]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chroms", default="1-22", help="Chromosome range/list, e.g. 1-22 or 1,2,22")
    p.add_argument("--skip-downloads", action="store_true", help="Skip GWAS raw data downloads")
    p.add_argument("--skip-pathways", action="store_true", help="Skip pathway aggregation stage")
    p.add_argument("--python", default=sys.executable, help="Python executable to use for sub-steps")
    return p.parse_args()


def run(cmd: list[str], desc: str) -> None:
    print(f"\n[{desc}] {' '.join(cmd)}", flush=True)
    env = {**__import__("os").environ, "PYTHONUNBUFFERED": "1"}
    r = subprocess.run(cmd, cwd=ROOT, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"FAILED ({r.returncode}): {desc}")


def concat_parquets(pattern: str, out_path: Path, required: bool = True) -> pd.DataFrame | None:
    files = sorted(out_path.parent.glob(pattern))
    if not files:
        if required:
            raise FileNotFoundError(f"No files matched {pattern}")
        return None
    frames = [pd.read_parquet(f) for f in files]
    merged = pd.concat(frames, ignore_index=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_parquet(out_path, index=False)
    print(f"Saved {out_path} ({len(merged):,} rows from {len(files)} chromosomes)")
    return merged


def assign_global_variant_split(master: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    uniq = master["variant_id"].dropna().unique()
    perm = rng.permutation(len(uniq))
    n_tr = int(0.70 * len(uniq))
    n_va = int(0.15 * len(uniq))
    assign: dict[str, str] = {}
    for i, idx in enumerate(perm):
        if i < n_tr:
            fold = "train"
        elif i < n_tr + n_va:
            fold = "val"
        else:
            fold = "test"
        assign[uniq[idx]] = fold
    master = master.copy()
    master["split_variant"] = master["variant_id"].map(assign)
    return master


def main() -> None:
    args = parse_args()
    chroms = parse_chrom_list(args.chroms)
    py = args.python

    print(f"Running genome-wide pipeline for chromosomes: {chroms}")

    # One-time prerequisite: harmonize PGS scores against all chr*.score.pvar files.
    # This overwrites the chr22-only harmonized TSVs so collect_variant_ids works
    # for every chromosome.
    run([py, "scripts/harmonize_pgs_genomewide.py"], "Harmonize PGS (genome-wide)")

    for chrom in chroms:
        if not args.skip_downloads:
            run([py, "scripts/download_panukbb_chrom.py", "--chrom", chrom], f"Pan-UKBB chr{chrom}")
            run([py, "scripts/download_bbj_chr22.py", "--chrom", chrom], f"BBJ chr{chrom}")
            run([py, "scripts/download_finngen_chr22.py", "--chrom", chrom], f"FinnGen chr{chrom}")

        run([py, "scripts/compute_af_features.py", "--chrom", chrom], f"AF features chr{chrom}")
        run([py, "scripts/compute_ld_features.py", "--chrom", chrom], f"LD features chr{chrom}")
        run([py, "scripts/compute_selection_features.py", "--chrom", chrom], f"Selection features chr{chrom}")
        run([py, "scripts/build_concordance_labels_multisource.py", "--chrom", chrom], f"Labels chr{chrom}")
        run([py, "scripts/build_master_table.py", "--chrom", chrom], f"Master table chr{chrom}")
        if not args.skip_pathways:
            run([py, "scripts/build_pathway_risk.py", "--chrom", chrom], f"Pathway risk chr{chrom}")
            run([py, "scripts/compute_constraint_features.py", "--chrom", chrom], f"Constraint features chr{chrom}")
            run([py, "scripts/build_master_table.py", "--chrom", chrom], f"Master table refresh chr{chrom}")

    labels_out = ROOT / "data/labels/gwas_concordance_labels_multisource_genomewide.parquet"
    concat_parquets("gwas_concordance_labels_multisource.chr*.parquet", labels_out)

    master_out = ROOT / "data/modeling/master_variant_table_genomewide.parquet"
    master = concat_parquets("master_variant_table.chr*.parquet", master_out)
    if master is not None:
        master = assign_global_variant_split(master)
        master.to_parquet(master_out, index=False)
        print(
            f"Reassigned global split_variant on genome-wide table "
            f"({master['variant_id'].nunique():,} variants)"
        )

    if not args.skip_pathways:
        concat_parquets(
            "pathway_risk_table.chr*.parquet",
            ROOT / "results/tables/pathway_risk_table_genomewide.parquet",
            required=False,
        )

    print("\nGenome-wide run complete.")


if __name__ == "__main__":
    main()
