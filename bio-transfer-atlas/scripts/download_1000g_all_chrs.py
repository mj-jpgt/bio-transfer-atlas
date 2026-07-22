"""
Download 1000 Genomes Phase 3 VCFs for all autosomes (chr1-22).
chr22 is skipped if already present.
Resumes partial downloads automatically.

Usage: python scripts/download_1000g_all_chrs.py
"""
import hashlib
import sys
import time
from pathlib import Path

import requests
from tqdm import tqdm

BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
TEMPLATE = "ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
CHRS = list(range(1, 23))

root = Path(__file__).resolve().parents[1]
vcf_dir = root / "data/raw/1000g/vcf"
vcf_dir.mkdir(parents=True, exist_ok=True)


def download_file(url: str, dest: Path, desc: str) -> bool:
    """Stream download with resume support. Returns True on success."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")

    # Resume if partial download exists
    resume_bytes = tmp.stat().st_size if tmp.exists() else 0
    headers = {"Range": f"bytes={resume_bytes}-"} if resume_bytes > 0 else {}

    try:
        r = requests.get(url, headers=headers, stream=True, timeout=60)
        if r.status_code == 416:  # Range not satisfiable = already complete
            tmp.rename(dest)
            return True
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0)) + resume_bytes

        mode = "ab" if resume_bytes > 0 else "wb"
        with open(tmp, mode) as f, tqdm(
            desc=desc,
            total=total,
            initial=resume_bytes,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            leave=False,
        ) as bar:
            for chunk in r.iter_content(chunk_size=1 << 20):  # 1 MB chunks
                f.write(chunk)
                bar.update(len(chunk))

        tmp.rename(dest)
        return True

    except Exception as e:
        print(f"\n  ERROR downloading {desc}: {e}", file=sys.stderr)
        return False


def main():
    print(f"Downloading 1000G Phase 3 autosomes chr1–22 → {vcf_dir}\n")

    failed = []
    for c in CHRS:
        fname = TEMPLATE.format(c=c)
        dest = vcf_dir / fname
        tbi_dest = vcf_dir / (fname + ".tbi")

        if dest.exists() and tbi_dest.exists():
            print(f"  chr{c:2d}: already present ({dest.stat().st_size / 1e6:.0f} MB) ✓")
            continue

        # VCF
        if not dest.exists():
            url = f"{BASE}/{fname}"
            ok = download_file(url, dest, f"chr{c} VCF")
            if not ok:
                failed.append(c)
                continue
        sz = dest.stat().st_size / 1e6
        print(f"  chr{c:2d}: downloaded {sz:.0f} MB ✓")

        # TBI index
        if not tbi_dest.exists():
            tbi_url = f"{BASE}/{fname}.tbi"
            download_file(tbi_url, tbi_dest, f"chr{c} .tbi")

        time.sleep(0.5)

    if failed:
        print(f"\nFailed chromosomes: {failed}")
        print("Re-run script to retry — partial files will be resumed.")
    else:
        print("\nAll chromosomes downloaded successfully.")


if __name__ == "__main__":
    main()
