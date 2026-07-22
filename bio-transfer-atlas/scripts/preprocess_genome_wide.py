"""
Genome-wide preprocessing pipeline for 1000G Phase 3 autosomes.
Runs after download_1000g_all_chrs.py completes.

Steps per chromosome:
  1. VCF → pgen  (biallelic SNPs only)
  2. QC pgen     (PCA dataset: geno 0.02, mind 0.02, MAF 0.01)
  3. Score pgen  (no MAF filter, for PGS scoring)
  4. Assign variant IDs  chr:pos:ref:alt

Genome-wide steps:
  5. LD prune each QC pgen, collect prune.in lists
  6. Merge all QC pgens into genome-wide pgen
  7. PCA on merged pruned variants (20 PCs)
  8. Save ancestry_pcs.parquet (replaces chr22-only version)

Usage: python scripts/preprocess_genome_wide.py [--chrs 1-22] [--threads 8]
"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parents[1]
plink2 = str(root / "tools/plink2/plink2.exe")
vcf_dir = root / "data/raw/1000g/vcf"
interim = root / "data/interim/1000g"
processed = root / "data/processed/1000g"
interim.mkdir(parents=True, exist_ok=True)
processed.mkdir(parents=True, exist_ok=True)

TEMPLATE = "ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"


def run(cmd: list, desc: str):
    print(f"  [{desc}]", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        # Print last 20 lines of stderr
        for line in r.stderr.strip().split("\n")[-20:]:
            print(f"    {line}")
        raise RuntimeError(f"FAILED: {desc}")
    return r


def vcf_to_pgen(chrom: int, threads: int):
    vcf = vcf_dir / TEMPLATE.format(c=chrom)
    out = interim / f"chr{chrom}"
    if (out.parent / f"chr{chrom}.pgen").exists():
        print(f"  chr{chrom}: pgen exists, skip")
        return
    run([plink2,
         "--vcf", str(vcf),
         "--max-alleles", "2",
         "--snps-only", "just-acgt",
         "--threads", str(threads),
         "--make-pgen",
         "--out", str(out)],
        f"chr{chrom} VCF→pgen")


def qc_pgen(chrom: int, threads: int):
    src = interim / f"chr{chrom}"
    out_pca = interim / f"chr{chrom}.qc.id"
    out_score = interim / f"chr{chrom}.score"

    # PCA dataset: MAF-filtered + unique IDs; --rm-dup removes same chr:pos:ref:alt duplicates
    if not (interim / f"chr{chrom}.qc.id.pgen").exists():
        run([plink2,
             "--pfile", str(src),
             "--geno", "0.02", "--mind", "0.02", "--maf", "0.01",
             "--set-all-var-ids", "@:#:$r:$a",
             "--new-id-max-allele-len", "1000",
             "--rm-dup", "exclude-mismatch",
             "--threads", str(threads),
             "--make-pgen",
             "--out", str(out_pca)],
            f"chr{chrom} QC+ID (PCA dataset)")

    # Scoring dataset: no MAF filter + unique IDs
    if not (interim / f"chr{chrom}.score.pgen").exists():
        run([plink2,
             "--pfile", str(src),
             "--geno", "0.02", "--mind", "0.02",
             "--set-all-var-ids", "@:#:$r:$a",
             "--new-id-max-allele-len", "1000",
             "--rm-dup", "exclude-mismatch",
             "--threads", str(threads),
             "--make-pgen",
             "--out", str(out_score)],
            f"chr{chrom} QC (scoring dataset, no MAF filter)")


def ld_prune(chrom: int, threads: int):
    src = interim / f"chr{chrom}.qc.id"
    out = interim / f"chr{chrom}.prune"
    if (interim / f"chr{chrom}.prune.prune.in").exists():
        return
    run([plink2,
         "--pfile", str(src),
         "--indep-pairwise", "200", "50", "0.2",
         "--threads", str(threads),
         "--out", str(out)],
        f"chr{chrom} LD prune")


def merge_and_pca(chrs: list, threads: int):
    """Merge all QC pgens, extract pruned variants, run PCA."""
    # Write pmerge list
    merge_list = interim / "pmerge_list.txt"
    lines = []
    for c in chrs:
        pfile = interim / f"chr{c}.qc.id"
        if (interim / f"chr{c}.qc.id.pgen").exists():
            lines.append(str(pfile))
    with open(merge_list, "w") as f:
        f.write("\n".join(lines[1:]))  # first file goes to --pfile, rest to --pmerge-list

    merged = interim / "genome_wide.qc"
    if not (interim / "genome_wide.qc.pgen").exists():
        run([plink2,
             "--pfile", lines[0],
             "--pmerge-list", str(merge_list),
             "--make-pgen",
             "--threads", str(threads),
             "--out", str(merged)],
            "Merge all chr pgens → genome_wide.qc")

    # Collect all prune.in lists
    prune_all = interim / "all_chrs.prune.in"
    if not prune_all.exists():
        with open(prune_all, "w") as out:
            for c in chrs:
                pin = interim / f"chr{c}.prune.prune.in"
                if pin.exists():
                    out.write(pin.read_text())

    # PCA
    pca_out = processed / "genome_wide_pca"
    if not (processed / "genome_wide_pca.eigenvec").exists():
        run([plink2,
             "--pfile", str(merged),
             "--extract", str(prune_all),
             "--pca", "20",
             "--threads", str(threads),
             "--out", str(pca_out)],
            "Genome-wide PCA (20 PCs)")

    print(f"\n  PCA outputs: {pca_out}.eigenvec / .eigenval")


def save_genome_wide_parquets():
    """Merge PCA eigenvec with population labels → parquet."""
    eigenvec = processed / "genome_wide_pca.eigenvec"
    if not eigenvec.exists():
        print("  PCA not yet complete, skipping parquet save.")
        return

    panel = pd.read_parquet(root / "data/processed/sample_metadata.parquet")
    pcs = pd.read_csv(eigenvec, sep="\t", header=0)
    pcs = pcs.rename(columns={"#IID": "sample_id"})
    if "FID" in pcs.columns:
        pcs = pcs.drop(columns=["FID"])
    pcs = pcs.merge(panel[["sample_id", "pop", "super_pop"]], on="sample_id", how="left")
    out = root / "data/processed/ancestry_pcs_genome_wide.parquet"
    pcs.to_parquet(out, index=False)
    print(f"  Saved: {out}  ({len(pcs)} samples × {len(pcs.columns)} cols)")

    eigenval = pd.read_csv(processed / "genome_wide_pca.eigenval", header=None, names=["eigenvalue"])
    eigenval["variance_explained"] = eigenval["eigenvalue"] / eigenval["eigenvalue"].sum()
    eigenval["PC"] = [f"PC{i+1}" for i in range(len(eigenval))]
    eigenval.to_parquet(processed / "genome_wide_pca_eigenvalues.parquet", index=False)
    print(f"  Top 5 PCs variance explained:")
    print(eigenval.head(5).to_string(index=False))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chrs", default="1-22",
                        help="Chromosomes to process, e.g. '1-22' or '1,2,3'")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip-download-check", action="store_true")
    args = parser.parse_args()

    # Parse chromosome list
    if "-" in args.chrs:
        start, end = args.chrs.split("-")
        chrs = list(range(int(start), int(end) + 1))
    else:
        chrs = [int(c) for c in args.chrs.split(",")]

    print(f"Processing chromosomes: {chrs}")
    print(f"Threads: {args.threads}\n")

    # Check all VCFs exist
    missing = [c for c in chrs if not (vcf_dir / TEMPLATE.format(c=c)).exists()]
    if missing and not args.skip_download_check:
        print(f"Missing VCFs for: chr{missing}")
        print("Run: python scripts/download_1000g_all_chrs.py  first")
        sys.exit(1)

    # Per-chromosome steps
    for c in chrs:
        print(f"\n=== chr{c} ===")
        vcf_to_pgen(c, args.threads)
        qc_pgen(c, args.threads)
        ld_prune(c, args.threads)

    # Genome-wide steps
    print("\n=== Genome-wide merge + PCA ===")
    merge_and_pca(chrs, args.threads)
    save_genome_wide_parquets()
    print("\nDone.")


if __name__ == "__main__":
    main()
