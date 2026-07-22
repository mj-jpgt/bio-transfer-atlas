"""
Conservative chr7+ score-pfile repair: PLINK work on local disk, low memory, optional qc.id fallback.

Use when OneDrive stubs or RAM pressure break the standard preprocess path.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from bta_plink import plink2_bin  # noqa: E402

PLINK2 = plink2_bin()
VCF_DIR = ROOT / "data/raw/1000g/vcf_grch38"
INTERIM = ROOT / "data/interim/1000g_grch38"
VCF_TEMPLATE = "1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Rebuild chr{N}.score.pgen conservatively.")
    p.add_argument("--chrom", required=True, help="Chromosome number, e.g. 7")
    p.add_argument("--memory-mb", type=int, default=2048, help="PLINK --memory (MB)")
    p.add_argument("--threads", type=int, default=1, help="PLINK --threads")
    p.add_argument("--min-free-gb", type=float, default=4.0, help="Abort if free RAM below this")
    p.add_argument(
        "--allow-qc-fallback",
        action="store_true",
        help="If VCF import fails, clone qc.id -> score",
    )
    p.add_argument(
        "--qc-fallback-only",
        action="store_true",
        help="Skip VCF import; clone qc.id -> score only (lowest RAM)",
    )
    p.add_argument("--work-root", default="", help="Local work dir (default: %TEMP%/bta_chr{N})")
    return p.parse_args()


def free_ram_gb() -> float:
    if sys.platform == "win32":
        try:
            import ctypes

            class MEMORYSTATUSEX(ctypes.Structure):
                _fields_ = [
                    ("dwLength", ctypes.c_ulong),
                    ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
                ]

            stat = MEMORYSTATUSEX()
            stat.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
            return stat.ullAvailPhys / (1024**3)
        except Exception:
            pass
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / (1024**2)
    except OSError:
        pass
    return 0.0


def plink(cmd: list[str], desc: str) -> None:
    print(f"  [plink] {desc}", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        err = (r.stderr or r.stdout or "").strip().split("\n")[-15:]
        for line in err:
            if line.strip():
                print(f"    {line}", flush=True)
        raise RuntimeError(f"PLINK failed: {desc}")


def score_pfile_ok(chrom: str, *, plink_check: bool = True) -> bool:
    pgen = INTERIM / f"chr{chrom}.score.pgen"
    pvar = INTERIM / f"chr{chrom}.score.pvar"
    psam = INTERIM / f"chr{chrom}.score.psam"
    if not all(p.exists() and p.stat().st_size > 1000 for p in (pgen, pvar, psam)):
        return False
    if not plink_check or free_ram_gb() < 1.5:
        return True
    out = INTERIM / f"_score_sanity_chr{chrom}"
    plink(
        [PLINK2, "--memory", "640", "--pfile", str(INTERIM / f"chr{chrom}.score"),
         "--freq", "--out", str(out)],
        f"sanity chr{chrom}.score",
    )
    return (Path(str(out) + ".afreq")).exists()


def copy_pfile_prefix(src: Path, dst: Path) -> None:
    for ext in (".pgen", ".pvar", ".psam", ".log"):
        s = Path(str(src) + ext)
        d = Path(str(dst) + ext)
        if s.exists():
            shutil.copy2(s, d)
            print(f"  copied -> {d.name}", flush=True)


def rebuild_from_vcf(chrom: str, work: Path, memory_mb: int, threads: int) -> None:
    vcf_name = VCF_TEMPLATE.format(chrom=chrom)
    vcf_src = VCF_DIR / vcf_name
    if not vcf_src.exists():
        raise FileNotFoundError(f"Missing VCF: {vcf_src}")

    work.mkdir(parents=True, exist_ok=True)
    vcf_local = work / vcf_name
    tbi_local = work / (vcf_name + ".tbi")
    if not vcf_local.exists():
        print(f"  Copying VCF to local work dir ({vcf_src.stat().st_size / 1e9:.2f} GB) ...", flush=True)
        shutil.copy2(vcf_src, vcf_local)
    tbi_src = Path(str(vcf_src) + ".tbi")
    if tbi_src.exists() and not tbi_local.exists():
        shutil.copy2(tbi_src, tbi_local)

    raw = work / f"chr{chrom}"
    score = work / f"chr{chrom}.score"

    if not Path(str(raw) + ".pgen").exists():
        plink(
            [
                PLINK2, "--memory", str(memory_mb),
                "--vcf", str(vcf_local),
                "--max-alleles", "2", "--snps-only", "just-acgt",
                "--threads", str(threads),
                "--make-pgen", "--out", str(raw),
            ],
            f"chr{chrom} VCF->pgen (local)",
        )

    if not Path(str(score) + ".pgen").exists():
        plink(
            [
                PLINK2, "--memory", str(memory_mb),
                "--pfile", str(raw),
                "--geno", "0.02", "--mind", "0.02",
                "--set-all-var-ids", "@:#:$r:$a",
                "--new-id-max-allele-len", "1000",
                "--rm-dup", "exclude-mismatch",
                "--threads", str(threads),
                "--make-pgen", "--out", str(score),
            ],
            f"chr{chrom} QC scoring pfile (local)",
        )

    print("  Copying score pfile to project interim ...", flush=True)
    copy_pfile_prefix(score, INTERIM / f"chr{chrom}.score")


def rebuild_from_qc_fallback(chrom: str, memory_mb: int, threads: int) -> None:
    qc = INTERIM / f"chr{chrom}.qc.id"
    if not Path(str(qc) + ".pgen").exists():
        raise FileNotFoundError(f"No qc.id pfile for chr{chrom}; cannot fall back")
    out = INTERIM / f"chr{chrom}.score"
    plink(
        [
            PLINK2, "--memory", str(memory_mb),
            "--pfile", str(qc),
            "--make-pgen",
            "--threads", str(threads),
            "--out", str(out),
        ],
        f"chr{chrom} qc.id -> score (fallback)",
    )
    print(
        "  WARNING: score pfile cloned from qc.id (MAF>=0.01 subset). "
        "Slightly fewer variants than full scoring build.",
        flush=True,
    )


def main() -> None:
    args = parse_args()
    chrom = str(int(args.chrom))
    free_gb = free_ram_gb()
    print(f"Free RAM: {free_gb:.1f} GB", flush=True)
    if free_gb < args.min_free_gb:
        raise SystemExit(
            f"Need >= {args.min_free_gb} GB free RAM (have {free_gb:.1f}). "
            "Close other apps and retry."
        )

    if score_pfile_ok(chrom, plink_check=free_ram_gb() >= 1.5):
        print(f"chr{chrom}.score.pgen OK — nothing to do.", flush=True)
        return

    if args.work_root:
        work = Path(args.work_root)
    else:
        tmp = os.environ.get("TMPDIR") or os.environ.get("TEMP") or "/tmp"
        work = Path(tmp) / f"bta_chr{chrom}"
    print(f"Work dir: {work}", flush=True)

    try:
        if args.qc_fallback_only:
            rebuild_from_qc_fallback(chrom, args.memory_mb, args.threads)
        else:
            rebuild_from_vcf(chrom, work, args.memory_mb, args.threads)
    except Exception as exc:
        print(f"VCF rebuild failed: {exc}", flush=True)
        if not args.allow_qc_fallback and not args.qc_fallback_only:
            raise
        rebuild_from_qc_fallback(chrom, args.memory_mb, args.threads)

    if not score_pfile_ok(chrom, plink_check=free_ram_gb() >= 1.5):
        raise RuntimeError(f"chr{chrom}.score still not readable after rebuild")
    print(f"chr{chrom}.score.pgen ready.", flush=True)


if __name__ == "__main__":
    main()
