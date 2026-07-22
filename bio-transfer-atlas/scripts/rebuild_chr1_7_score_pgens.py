"""

Phase A5: Orchestrate chr1-7 score.pgen rebuild (VCF or qc-fallback).

Supports --jobs N concurrent chromosome rebuilds (Lambda profile).

"""

from __future__ import annotations



import argparse

import os

import subprocess

import sys

from concurrent.futures import ThreadPoolExecutor, as_completed

from pathlib import Path



ROOT = Path(__file__).resolve().parents[1]

INTERIM = ROOT / "data/interim/1000g_grch38"

VCF_DIR = ROOT / "data/raw/1000g/vcf_grch38"

PY = sys.executable





def parse_args() -> argparse.Namespace:

    p = argparse.ArgumentParser()

    p.add_argument("--chroms", default="1,2,3,4,5,6,7")

    p.add_argument("--memory-mb", type=int, default=int(os.environ.get("BTA_PLINK_MEMORY_MB", "640")))

    p.add_argument("--threads", type=int, default=int(os.environ.get("BTA_PLINK_THREADS", "1")))

    p.add_argument("--jobs", type=int, default=int(os.environ.get("BTA_REBUILD_JOBS", "1")))

    p.add_argument("--qc-fallback-only", action="store_true")

    return p.parse_args()





def rebuild_one(chrom: str, memory_mb: int, threads: int, qc_fallback_only: bool) -> dict:

    score = INTERIM / f"chr{chrom}.score.pgen"

    if score.exists() and score.stat().st_size > 1000:

        print(f"chr{chrom}: score.pgen already present", flush=True)

        return {"chrom": chrom, "status": "exists"}



    qc = INTERIM / f"chr{chrom}.qc.pgen"

    vcf_long = VCF_DIR / (

        f"1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

    )

    vcf_short = VCF_DIR / f"chr{chrom}.vcf.gz"

    cmd = [

        PY,

        "scripts/rebuild_score_pfile_conservative.py",

        "--chrom",

        chrom,

        "--memory-mb",

        str(memory_mb),

        "--threads",

        str(threads),

        "--min-free-gb",

        "0",

        "--work-root",

        str(Path(os.environ.get("TMPDIR", "/tmp")) / f"bta_chr{chrom}"),

    ]

    if qc_fallback_only or (qc.exists() and not vcf_long.exists() and not vcf_short.exists()):

        if not qc.exists():

            print(f"chr{chrom}: MISSING qc.pgen and VCF — cannot rebuild", flush=True)

            return {"chrom": chrom, "status": "blocked_no_vcf_no_qc"}

        cmd.append("--qc-fallback-only")

        print(f"chr{chrom}: qc-fallback rebuild", flush=True)

    else:

        print(f"chr{chrom}: attempting VCF rebuild", flush=True)

        cmd.append("--allow-qc-fallback")



    r = subprocess.run(cmd, cwd=ROOT)

    status = "ok" if r.returncode == 0 and score.exists() else "failed"

    print(f"chr{chrom}: {status}", flush=True)

    return {"chrom": chrom, "status": status}





def main() -> None:

    args = parse_args()

    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]

    INTERIM.mkdir(parents=True, exist_ok=True)

    jobs = max(1, args.jobs)

    report: list[dict] = []



    if jobs == 1:

        for chrom in chroms:

            report.append(rebuild_one(chrom, args.memory_mb, args.threads, args.qc_fallback_only))

    else:

        with ThreadPoolExecutor(max_workers=jobs) as ex:

            futs = {

                ex.submit(rebuild_one, c, args.memory_mb, args.threads, args.qc_fallback_only): c

                for c in chroms

            }

            for fut in as_completed(futs):

                report.append(fut.result())



    # stable chrom order

    order = {c: i for i, c in enumerate(chroms)}

    report.sort(key=lambda r: order.get(str(r["chrom"]), 99))



    out = ROOT / "results/tables/score_pgen_chr1_7_rebuild_status.csv"

    out.parent.mkdir(parents=True, exist_ok=True)

    import pandas as pd



    pd.DataFrame(report).to_csv(out, index=False)

    print(f"Saved {out}")

    blocked = [r for r in report if str(r["status"]).startswith("blocked") or r["status"] == "failed"]

    if blocked:

        print("NOTE: some chroms blocked — restore VCFs from IGSR/HF then re-run")





if __name__ == "__main__":

    main()


