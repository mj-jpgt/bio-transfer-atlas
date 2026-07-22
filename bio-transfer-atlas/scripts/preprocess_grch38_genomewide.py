"""
Preprocess 1000G GRCh38 chr1-22 into per-chromosome scoring pfiles used by atlas features.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PLINK2 = str(ROOT / "tools/plink2/plink2.exe")
PLINK_MEMORY_MB = 4096
VCF_DIR = ROOT / "data/raw/1000g/vcf_grch38"
INTERIM = ROOT / "data/interim/1000g_grch38"
PROCESSED = ROOT / "data/processed/1000g_grch38"
INTERIM.mkdir(parents=True, exist_ok=True)
PROCESSED.mkdir(parents=True, exist_ok=True)


def parse_chrom_list(spec: str) -> list[int]:
    s = spec.strip()
    if "-" in s:
        a, b = s.split("-", 1)
        return list(range(int(a), int(b) + 1))
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrs", default="1-22", help="Chromosome range/list, e.g. 1-22 or 1,2,22")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--skip-merge", action="store_true", help="Skip genome-wide merge/PCA (single-chr rebuild)")
    p.add_argument(
        "--vcf-template",
        default="1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
        help="Per-chromosome VCF filename template",
    )
    return p.parse_args()


def run(cmd: list[str], desc: str) -> None:
    print(f"  [{desc}]", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        for line in r.stderr.strip().split("\n")[-20:]:
            if line.strip():
                print(f"    {line}")
        raise RuntimeError(f"FAILED: {desc}")


def vcf_to_pgen(chrom: int, vcf_template: str, threads: int) -> None:
    out = INTERIM / f"chr{chrom}"
    if (INTERIM / f"chr{chrom}.pgen").exists():
        print(f"  chr{chrom}: raw pgen exists, skip")
        return
    vcf = VCF_DIR / vcf_template.format(chrom=chrom)
    if not vcf.exists():
        raise FileNotFoundError(f"Missing VCF for chr{chrom}: {vcf}")
    run(
        [
            PLINK2,
            "--memory",
            str(PLINK_MEMORY_MB),
            "--vcf",
            str(vcf),
            "--max-alleles",
            "2",
            "--snps-only",
            "just-acgt",
            "--threads",
            str(threads),
            "--make-pgen",
            "--out",
            str(out),
        ],
        f"chr{chrom} VCF->pgen",
    )


def qc_pgens(chrom: int, threads: int) -> None:
    src = INTERIM / f"chr{chrom}"
    out_qc = INTERIM / f"chr{chrom}.qc.id"
    out_score = INTERIM / f"chr{chrom}.score"
    if not (INTERIM / f"chr{chrom}.qc.id.pgen").exists():
        run(
            [
                PLINK2,
                "--memory",
                str(PLINK_MEMORY_MB),
                "--pfile",
                str(src),
                "--geno",
                "0.02",
                "--mind",
                "0.02",
                "--maf",
                "0.01",
                "--set-all-var-ids",
                "@:#:$r:$a",
                "--new-id-max-allele-len",
                "1000",
                "--rm-dup",
                "exclude-mismatch",
                "--threads",
                str(threads),
                "--make-pgen",
                "--out",
                str(out_qc),
            ],
            f"chr{chrom} QC+ID (PCA dataset)",
        )
    if not (INTERIM / f"chr{chrom}.score.pgen").exists():
        run(
            [
                PLINK2,
                "--memory",
                str(PLINK_MEMORY_MB),
                "--pfile",
                str(src),
                "--geno",
                "0.02",
                "--mind",
                "0.02",
                "--set-all-var-ids",
                "@:#:$r:$a",
                "--new-id-max-allele-len",
                "1000",
                "--rm-dup",
                "exclude-mismatch",
                "--threads",
                str(threads),
                "--make-pgen",
                "--out",
                str(out_score),
            ],
            f"chr{chrom} QC (scoring dataset, no MAF filter)",
        )


def ld_prune(chrom: int, threads: int) -> None:
    src = INTERIM / f"chr{chrom}.qc.id"
    out = INTERIM / f"chr{chrom}.prune"
    if (INTERIM / f"chr{chrom}.prune.prune.in").exists():
        return
    run(
        [
            PLINK2,
            "--memory",
            str(PLINK_MEMORY_MB),
            "--pfile",
            str(src),
            "--indep-pairwise",
            "200",
            "50",
            "0.2",
            "--threads",
            str(threads),
            "--out",
            str(out),
        ],
        f"chr{chrom} LD prune",
    )


def merge_and_pca(chrs: list[int], threads: int) -> None:
    merge_list = INTERIM / "pmerge_list_grch38.txt"
    lines = []
    for c in chrs:
        pfile = INTERIM / f"chr{c}.qc.id"
        if (INTERIM / f"chr{c}.qc.id.pgen").exists():
            lines.append(str(pfile))
    if not lines:
        raise RuntimeError("No chr*.qc.id pgen files found for merge")
    merge_list.write_text("\n".join(lines[1:]), encoding="utf-8")

    merged = INTERIM / "genome_wide_grch38.qc"
    if not (INTERIM / "genome_wide_grch38.qc.pgen").exists():
        run(
            [
                PLINK2,
                "--memory",
                str(PLINK_MEMORY_MB),
                "--pfile",
                lines[0],
                "--pmerge-list",
                str(merge_list),
                "--make-pgen",
                "--threads",
                str(threads),
                "--out",
                str(merged),
            ],
            "Merge all chr QC pgens (GRCh38)",
        )

    prune_all = INTERIM / "all_chrs_grch38.prune.in"
    if not prune_all.exists():
        with open(prune_all, "w", encoding="utf-8") as out:
            for c in chrs:
                pin = INTERIM / f"chr{c}.prune.prune.in"
                if pin.exists():
                    out.write(pin.read_text(encoding="utf-8"))

    pca_out = PROCESSED / "genome_wide_grch38_pca"
    if not (PROCESSED / "genome_wide_grch38_pca.eigenvec").exists():
        run(
            [
                PLINK2,
                "--memory",
                str(PLINK_MEMORY_MB),
                "--pfile",
                str(merged),
                "--extract",
                str(prune_all),
                "--pca",
                "20",
                "--threads",
                str(threads),
                "--out",
                str(pca_out),
            ],
            "Genome-wide PCA GRCh38 (20 PCs)",
        )


def save_parquets() -> None:
    eigenvec = PROCESSED / "genome_wide_grch38_pca.eigenvec"
    if not eigenvec.exists():
        print("  PCA not found; skipping parquet export.")
        return
    panel = pd.read_parquet(ROOT / "data/processed/sample_metadata_grch38.parquet")
    pcs = pd.read_csv(eigenvec, sep="\t", header=0)
    pcs = pcs.rename(columns={"#IID": "sample_id"})
    if "FID" in pcs.columns:
        pcs = pcs.drop(columns=["FID"])
    pcs = pcs.merge(
        panel[["sample_id", "pop", "super_pop"]],
        on="sample_id",
        how="left",
        validate="many_to_one",
    )
    out = ROOT / "data/processed/ancestry_pcs_grch38_genomewide.parquet"
    pcs.to_parquet(out, index=False)
    print(f"  Saved: {out} ({len(pcs):,} samples × {len(pcs.columns)} cols)")


def cleanup_raw(chrom: int, vcf_template: str) -> None:
    """Delete the raw pgen/pvar and source VCF once score+qc pgens are ready.

    This keeps disk usage low enough to process all 22 chromosomes without
    running out of space (~the raw pgen+pvar alone are 6-10 GB per chromosome).
    """
    # Only clean up if both downstream outputs already exist
    score_ok = (INTERIM / f"chr{chrom}.score.pgen").exists()
    qc_ok = (INTERIM / f"chr{chrom}.qc.id.pgen").exists()
    prune_ok = (INTERIM / f"chr{chrom}.prune.prune.in").exists()
    if not (score_ok and qc_ok and prune_ok):
        return

    for ext in (".pgen", ".pvar", ".psam", ".log"):
        raw = INTERIM / f"chr{chrom}{ext}"
        if raw.exists():
            raw.unlink()
            print(f"  cleaned: {raw.name}")

    vcf = VCF_DIR / vcf_template.format(chrom=chrom)
    for p in (vcf, Path(str(vcf) + ".tbi")):
        if p.exists():
            p.unlink()
            print(f"  cleaned: {p.name}")


def main() -> None:
    args = parse_args()
    chrs = parse_chrom_list(args.chrs)
    print(f"Processing GRCh38 chromosomes: {chrs}")
    print(f"VCF dir: {VCF_DIR}")

    for c in chrs:
        print(f"\n=== chr{c} ===")
        vcf_to_pgen(c, args.vcf_template, args.threads)
        qc_pgens(c, args.threads)
        ld_prune(c, args.threads)
        # Remove raw pgen + VCF immediately to keep disk usage bounded
        cleanup_raw(c, args.vcf_template)

    print("\n=== Genome-wide merge + PCA ===")
    if args.skip_merge:
        print("  skipped (--skip-merge)")
    else:
        merge_and_pca(chrs, args.threads)
        save_parquets()
    print("\nDone.")


if __name__ == "__main__":
    main()
