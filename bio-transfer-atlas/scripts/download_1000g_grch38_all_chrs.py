"""
Download 1000G high-coverage GRCh38 phased panel chr1-22 VCFs + TBI indexes.
Supports --jobs for parallel chromosome downloads (Lambda).
"""
from __future__ import annotations

import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
)
TEMPLATE = "1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data/raw/1000g/vcf_grch38"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_chrom_list(spec: str) -> list[str]:
    s = spec.strip()
    if "-" in s and "," not in s:
        a, b = s.split("-", 1)
        return [str(x) for x in range(int(a), int(b) + 1)]
    return [str(int(x.strip())) for x in s.split(",") if x.strip()]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chroms", default="1-22", help="Chromosome range/list, e.g. 1-22 or 1,2,22")
    p.add_argument(
        "--jobs",
        type=int,
        default=int(os.environ.get("BTA_DOWNLOAD_JOBS", "1")),
        help="Parallel chromosome downloads",
    )
    return p.parse_args()


def download_one(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"cached -> {dest.name}", flush=True)
        return
    tmp = Path(str(dest) + ".partial")
    start = tmp.stat().st_size if tmp.exists() else 0
    range_arg = f"{start}-" if start > 0 else None
    from bta_curl import curl_bin

    cmd = [curl_bin(), "-fLsS", "--retry", "3", "--retry-delay", "5"]
    if range_arg is not None:
        cmd.extend(["-C", "-", "-o", str(tmp), url])  # continue
    else:
        cmd.extend(["-o", str(tmp), url])
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Download failed for {url}: {(r.stderr or r.stdout)[-400:]}")
    tmp.rename(dest)
    print(f"downloaded -> {dest.name} ({dest.stat().st_size / 1e9:.2f} GB)", flush=True)


def download_chrom(chrom: str) -> str:
    vcf_name = TEMPLATE.format(chrom=chrom)
    vcf_url = f"{BASE}/{vcf_name}"
    tbi_url = f"{vcf_url}.tbi"
    download_one(vcf_url, OUT_DIR / vcf_name)
    download_one(tbi_url, OUT_DIR / f"{vcf_name}.tbi")
    return chrom


def main() -> None:
    args = parse_args()
    chroms = parse_chrom_list(args.chroms)
    jobs = max(1, args.jobs)
    print(f"Downloading 1000G GRCh38 phased VCFs for chromosomes: {chroms} (jobs={jobs})")
    if jobs == 1:
        for chrom in chroms:
            download_chrom(chrom)
    else:
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = [ex.submit(download_chrom, c) for c in chroms]
            for fut in as_completed(futs):
                print(f"chrom {fut.result()} complete", flush=True)
    print("Done.")


if __name__ == "__main__":
    main()
