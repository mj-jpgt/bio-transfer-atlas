"""
Phase 18.4: Re-score 1000G genotypes with intervention-modified PGS weights.

Parallel chromosome jobs + resume via existing .sscore files (PLINK is CPU/RAM;
GPU does not accelerate --score).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from intervention_common import (
    GW_INT_ROOT,
    GW_TAG,
    INTERVENTION_MODES,
    PGS_IDS,
    ROOT,
    available_score_chroms,
    pfile_for_chrom,
    plink2_bin,
)

SCORE_OUT = ROOT / "data/processed/scores_grch38_intervention_genomewide"
PANEL_PATH = ROOT / "data/processed/sample_metadata_grch38.parquet"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="PLINK2 --score per chr, sum per intervention mode.")
    p.add_argument("--intervention-root", default=str(GW_INT_ROOT))
    p.add_argument("--out-dir", default=str(SCORE_OUT))
    p.add_argument("--chroms", default="", help="Comma list; default available score pgens")
    p.add_argument("--modes", default=",".join(INTERVENTION_MODES))
    p.add_argument("--pgs-ids", default=",".join(PGS_IDS))
    p.add_argument("--memory-mb", type=int, default=8192)
    p.add_argument("--threads", type=int, default=4, help="PLINK threads per job")
    p.add_argument("--jobs", type=int, default=8, help="Parallel score jobs")
    p.add_argument("--tag", default=GW_TAG, help="Suffix for score_matrix_{mode}_{tag}.parquet")
    return p.parse_args()


def parse_chrom_list(spec: str) -> list[str]:
    if not spec.strip():
        return available_score_chroms()
    return [str(int(x.strip())) for x in spec.split(",") if x.strip()]


def run_plink_score(
    pfile: Path,
    score_tsv: Path,
    out_prefix: Path,
    memory_mb: int,
    threads: int,
) -> Path | None:
    sscore = Path(str(out_prefix) + ".sscore")
    if sscore.exists() and sscore.stat().st_size > 0:
        return sscore
    cmd = [
        plink2_bin(),
        "--threads",
        str(max(1, threads)),
        "--pfile",
        str(pfile),
        "--memory",
        str(memory_mb),
        "--score",
        str(score_tsv),
        "1",
        "2",
        "3",
        "header",
        "cols=+scoresums",
        "--out",
        str(out_prefix),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  SCORE FAILED {score_tsv.name} {out_prefix.name}: {r.stderr[-300:]}", flush=True)
        return None
    return sscore if sscore.exists() else None


def sum_chr_sscores(sscore_paths: list[Path]) -> pd.Series:
    total: pd.Series | None = None
    for path in sscore_paths:
        df = pd.read_csv(path, sep="\t", usecols=["#IID", "SCORE1_SUM"])
        s = df.set_index("#IID")["SCORE1_SUM"]
        total = s if total is None else total.add(s, fill_value=0.0)
    return total if total is not None else pd.Series(dtype=float)


def main() -> None:
    args = parse_args()
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    pgs_ids = [p.strip() for p in args.pgs_ids.split(",") if p.strip()]
    chroms = parse_chrom_list(args.chroms)
    int_root = Path(args.intervention_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    panel = pd.read_parquet(PANEL_PATH)
    work = out_dir / "_sscore_work"
    work.mkdir(parents=True, exist_ok=True)
    jobs = max(1, args.jobs)

    if not chroms:
        raise SystemExit("No score pgens available")

    print(
        f"Chromosomes: {chroms} jobs={jobs} threads={args.threads} mem={args.memory_mb}",
        flush=True,
    )

    # Flatten all scoring tasks for max CPU/RAM utilization
    tasks: list[tuple[str, str, str, Path, Path]] = []
    for mode in modes:
        for pgs_id in pgs_ids:
            score_tsv = int_root / pgs_id / f"{mode}.tsv"
            if not score_tsv.exists():
                print(f"  skip {pgs_id}/{mode}: missing {score_tsv}", flush=True)
                continue
            for chrom in chroms:
                pfile = pfile_for_chrom(chrom)
                if not Path(str(pfile) + ".pgen").exists():
                    continue
                out_prefix = work / f"{mode}_{pgs_id}_chr{chrom}"
                tasks.append((mode, pgs_id, chrom, score_tsv, out_prefix))

    print(f"Total PLINK score tasks: {len(tasks)}", flush=True)
    done = 0
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {
            ex.submit(
                run_plink_score,
                pfile_for_chrom(chrom),
                score_tsv,
                out_prefix,
                args.memory_mb,
                args.threads,
            ): (mode, pgs_id, chrom)
            for mode, pgs_id, chrom, score_tsv, out_prefix in tasks
        }
        for fut in as_completed(futs):
            mode, pgs_id, chrom = futs[fut]
            sscore = fut.result()
            done += 1
            if sscore and done % 25 == 0:
                print(f"  progress {done}/{len(tasks)} (last {mode}/{pgs_id}/chr{chrom})", flush=True)

    print("Aggregating score matrices ...", flush=True)
    for mode in modes:
        pgs_series: dict[str, pd.Series] = {}
        for pgs_id in pgs_ids:
            chr_paths = []
            for chrom in chroms:
                sscore = work / f"{mode}_{pgs_id}_chr{chrom}.sscore"
                if sscore.exists() and sscore.stat().st_size > 0:
                    chr_paths.append(sscore)
            if not chr_paths:
                continue
            pgs_series[pgs_id] = sum_chr_sscores(chr_paths)
            print(f"  {mode} {pgs_id}: {len(chr_paths)} chroms", flush=True)

        if not pgs_series:
            print(f"  No scores for mode {mode}", flush=True)
            continue

        merged = None
        for pgs_id, series in pgs_series.items():
            col_df = series.rename(pgs_id).reset_index().rename(columns={"#IID": "sample_id"})
            merged = col_df if merged is None else merged.merge(col_df, on="sample_id", how="outer")

        merged = merged.merge(panel, on="sample_id", how="left")
        out_path = out_dir / f"score_matrix_{mode}_{args.tag}.parquet"
        merged.to_parquet(out_path, index=False)
        score_cols = [c for c in merged.columns if c.startswith("PGS")]
        print(f"  -> {out_path} ({len(merged)} samples x {len(score_cols)} PGS)", flush=True)

    print("SCORE_INTERVENTION_DONE", flush=True)


if __name__ == "__main__":
    main()
