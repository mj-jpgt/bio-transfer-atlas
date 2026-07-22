"""
Genome-wide PGS scoring: PLINK2 --score per chromosome, sum SCORE1_SUM additively.

Writes data/processed/scores_grch38/score_matrix_grch38_genomewide.parquet
Supports parallel chrom jobs and skips existing .sscore files (resume-friendly).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from intervention_common import PGS_IDS, PLINK2, ROOT, available_score_chroms, harmonized_path

INTERIM = ROOT / "data/interim/1000g_grch38"
SCORE_OUT = ROOT / "data/processed/scores_grch38"
PANEL_PATH = ROOT / "data/processed/sample_metadata_grch38.parquet"
DEFAULT_OUT = SCORE_OUT / "score_matrix_grch38_genomewide.parquet"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Sum per-chromosome PLINK PGS scores genome-wide.")
    p.add_argument("--chroms", default="", help="Comma list or 1-22; default all with score.pgen")
    p.add_argument("--pgs-ids", default=",".join(PGS_IDS))
    p.add_argument("--memory-mb", type=int, default=int(os.environ.get("BTA_PLINK_MEMORY_MB", "2048")))
    p.add_argument("--threads", type=int, default=int(os.environ.get("BTA_PLINK_THREADS", "8")))
    p.add_argument(
        "--jobs",
        type=int,
        default=int(os.environ.get("BTA_SCORE_JOBS", "4")),
        help="Parallel chromosome score jobs",
    )
    p.add_argument("--out", default=str(DEFAULT_OUT))
    p.add_argument("--work-dir", default="", help="Temp dir for .sscore files")
    return p.parse_args()


def parse_chrom_list(spec: str) -> list[str]:
    if not spec.strip():
        return available_score_chroms()
    s = spec.strip()
    if "-" in s and "," not in s:
        a, b = s.split("-", 1)
        return [str(x) for x in range(int(a), int(b) + 1)]
    return [str(int(x.strip())) for x in s.split(",") if x.strip()]


def run_plink_score(
    pfile: Path, score_tsv: Path, out_prefix: Path, memory_mb: int, threads: int
) -> Path | None:
    sscore = Path(str(out_prefix) + ".sscore")
    if sscore.exists() and sscore.stat().st_size > 100:
        return sscore
    cmd = [
        PLINK2,
        "--pfile",
        str(pfile),
        "--memory",
        str(memory_mb),
        "--threads",
        str(threads),
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
        err = (r.stderr or r.stdout or "").strip().split("\n")[-5:]
        print(f"  FAILED {out_prefix.name}: {' | '.join(err)}", flush=True)
        return None
    return sscore if sscore.exists() else None


def sum_chr_scores(sscore_paths: list[Path]) -> pd.Series:
    total: pd.Series | None = None
    for path in sscore_paths:
        df = pd.read_csv(path, sep="\t", usecols=["#IID", "SCORE1_SUM"])
        s = df.set_index("#IID")["SCORE1_SUM"]
        total = s if total is None else total.add(s, fill_value=0.0)
    if total is None:
        return pd.Series(dtype=float)
    return total


def main() -> None:
    args = parse_args()
    chroms = parse_chrom_list(args.chroms)
    pgs_ids = [p.strip() for p in args.pgs_ids.split(",") if p.strip()]
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    work = Path(args.work_dir) if args.work_dir else SCORE_OUT / "_sscore_work"
    work.mkdir(parents=True, exist_ok=True)
    jobs = max(1, args.jobs)

    if not chroms:
        raise SystemExit("No score pgens found under data/interim/1000g_grch38/")

    print(
        f"Chromosomes: {chroms} jobs={jobs} threads={args.threads} mem={args.memory_mb}",
        flush=True,
    )
    panel = pd.read_parquet(PANEL_PATH)
    merged: pd.DataFrame | None = None

    for pgs_id in pgs_ids:
        score_tsv = harmonized_path(pgs_id)
        if not score_tsv.exists():
            print(f"skip {pgs_id}: missing {score_tsv}", flush=True)
            continue

        tasks = []
        for chrom in chroms:
            pfile = INTERIM / f"chr{chrom}.score"
            if not Path(str(pfile) + ".pgen").exists():
                print(f"  skip chr{chrom}: no score pgen", flush=True)
                continue
            out_prefix = work / f"{pgs_id}.chr{chrom}"
            tasks.append((chrom, pfile, out_prefix))

        chr_sscores: list[Path] = []
        if jobs == 1:
            for chrom, pfile, out_prefix in tasks:
                sscore = run_plink_score(
                    pfile, score_tsv, out_prefix, args.memory_mb, args.threads
                )
                if sscore:
                    chr_sscores.append(sscore)
                    print(f"  {pgs_id} chr{chrom} scored", flush=True)
        else:
            with ThreadPoolExecutor(max_workers=jobs) as ex:
                futs = {
                    ex.submit(
                        run_plink_score,
                        pfile,
                        score_tsv,
                        out_prefix,
                        args.memory_mb,
                        args.threads,
                    ): chrom
                    for chrom, pfile, out_prefix in tasks
                }
                for fut in as_completed(futs):
                    chrom = futs[fut]
                    sscore = fut.result()
                    if sscore:
                        chr_sscores.append(sscore)
                        print(f"  {pgs_id} chr{chrom} scored", flush=True)

        if not chr_sscores:
            print(f"  no scores for {pgs_id}", flush=True)
            continue

        summed = sum_chr_scores(chr_sscores)
        col_df = summed.rename(pgs_id).reset_index().rename(columns={"#IID": "sample_id"})
        merged = col_df if merged is None else merged.merge(col_df, on="sample_id", how="outer")

    if merged is None:
        raise SystemExit("No PGS scores produced")

    merged = merged.merge(panel, on="sample_id", how="left")
    # Preserve previously scored PGS columns when expanding the matrix
    if out_path.exists():
        try:
            old = pd.read_parquet(out_path)
            keep = [c for c in old.columns if c.startswith("PGS") and c not in merged.columns]
            if keep:
                merged = merged.merge(old[["sample_id"] + keep], on="sample_id", how="left")
                print(f"  merged {len(keep)} existing PGS columns from prior matrix", flush=True)
        except Exception as exc:
            print(f"  warn: could not merge prior matrix: {exc}", flush=True)
    merged.to_parquet(out_path, index=False)
    score_cols = [c for c in merged.columns if c.startswith("PGS")]
    print(f"\nSaved {out_path} ({len(merged)} samples x {len(score_cols)} PGS)", flush=True)


if __name__ == "__main__":
    main()
