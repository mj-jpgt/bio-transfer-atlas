"""
FAIRGEN-Open Stage 2: GRCh38 chr22 Preprocessing (Smoke Test)
==============================================================
Input:  1kGP_high_coverage_Illumina chr22 phased panel (GRCh38, 3202 samples)
Output:
  data/interim/1000g_grch38/chr22.pgen          - raw converted
  data/interim/1000g_grch38/chr22.qc.pgen       - PCA dataset (MAF>1%)
  data/interim/1000g_grch38/chr22.score.pgen    - scoring dataset (no MAF filter)
  data/processed/1000g_grch38/chr22_pca.*       - PCA eigenvec/eigenval
  data/processed/sample_metadata_grch38.parquet - 3202 samples + pop labels
  data/processed/ancestry_pcs_chr22.parquet     - chr22-only PCA (smoke test)

⚠️  CHR22 SMOKE TEST — not genome-wide. PCA variance inflated.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parents[1]
plink2 = str(root / "tools/plink2/plink2.exe")

vcf_dir   = root / "data/raw/1000g/vcf_grch38"
meta_dir  = root / "data/raw/1000g/metadata"
interim   = root / "data/interim/1000g_grch38"
processed = root / "data/processed/1000g_grch38"
interim.mkdir(parents=True, exist_ok=True)
processed.mkdir(parents=True, exist_ok=True)

VCF = vcf_dir / "chr22.vcf.gz"
THREADS = 8


def run(cmd: list, desc: str):
    print(f"  [{desc}]", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        for line in r.stderr.strip().split("\n")[-20:]:
            if line.strip():
                print(f"    {line}")
        raise RuntimeError(f"FAILED: {desc}")
    return r


def step1_vcf_to_pgen():
    out = interim / "chr22"
    if (interim / "chr22.pgen").exists():
        print("  chr22.pgen exists, skip")
        return
    run([plink2,
         "--vcf", str(VCF),
         "--max-alleles", "2",
         "--snps-only", "just-acgt",
         "--threads", str(THREADS),
         "--make-pgen",
         "--out", str(out)],
        "chr22 VCF→pgen")
    # Report variant/sample counts
    pvar = pd.read_csv(str(out) + ".pvar", sep="\t", comment="#", header=None,
                       usecols=[0,1], nrows=5)
    psam = pd.read_csv(str(out) + ".psam", sep="\t", nrows=3)
    n_vars = sum(1 for _ in open(str(out) + ".pvar")) - 1
    n_samp = sum(1 for _ in open(str(out) + ".psam")) - 1
    print(f"    variants: {n_vars:,}   samples: {n_samp:,}")


def step2_qc_pca_pgen():
    out = interim / "chr22.qc"
    if (interim / "chr22.qc.pgen").exists():
        print("  chr22.qc.pgen exists, skip")
        return
    run([plink2,
         "--pfile", str(interim / "chr22"),
         "--geno", "0.02", "--mind", "0.02", "--maf", "0.01",
         "--set-all-var-ids", "@:#:$r:$a",
         "--new-id-max-allele-len", "1000",
         "--rm-dup", "exclude-mismatch",
         "--threads", str(THREADS),
         "--make-pgen",
         "--out", str(out)],
        "chr22 QC+ID (PCA dataset, MAF>1%)")


def step3_score_pgen():
    out = interim / "chr22.score"
    if (interim / "chr22.score.pgen").exists():
        print("  chr22.score.pgen exists, skip")
        return
    run([plink2,
         "--pfile", str(interim / "chr22"),
         "--geno", "0.02", "--mind", "0.02",
         "--set-all-var-ids", "@:#:$r:$a",
         "--new-id-max-allele-len", "1000",
         "--rm-dup", "exclude-mismatch",
         "--threads", str(THREADS),
         "--make-pgen",
         "--out", str(out)],
        "chr22 QC (scoring dataset, no MAF filter)")


def step4_ld_prune():
    prune_out = interim / "chr22.prune"
    if (interim / "chr22.prune.prune.in").exists():
        print("  prune.in exists, skip")
        return
    run([plink2,
         "--pfile", str(interim / "chr22.qc"),
         "--indep-pairwise", "200", "50", "0.2",
         "--threads", str(THREADS),
         "--out", str(prune_out)],
        "chr22 LD prune")


def step5_pca():
    pca_out = processed / "chr22_pca"
    if (processed / "chr22_pca.eigenvec").exists():
        print("  PCA eigenvec exists, skip")
        return
    run([plink2,
         "--pfile", str(interim / "chr22.qc"),
         "--extract", str(interim / "chr22.prune.prune.in"),
         "--pca", "20",
         "--threads", str(THREADS),
         "--out", str(pca_out)],
        "chr22 PCA (20 PCs)")


def step6_save_parquets():
    """Merge PCA + sample panel → parquet files."""
    panel_path = meta_dir / "integrated_call_samples_v3.20130502.ALL.panel"
    eigenvec   = processed / "chr22_pca.eigenvec"
    eigenval   = processed / "chr22_pca.eigenval"

    panel = pd.read_csv(panel_path, sep="\t")
    panel.columns = panel.columns.str.strip()

    # Normalise column names — Phase3 panel uses 'sample', 'pop', 'super_pop', 'gender'
    rename = {}
    for c in panel.columns:
        lc = c.lower()
        if lc in ("sample", "individual_id", "#sample"): rename[c] = "sample_id"
        elif lc == "pop": rename[c] = "pop"
        elif lc in ("super_pop", "super population code"): rename[c] = "super_pop"
        elif lc in ("gender", "sex"): rename[c] = "gender"
    panel = panel.rename(columns=rename)
    panel = panel[["sample_id", "pop", "super_pop", "gender"]].drop_duplicates()
    out_meta = root / "data/processed/sample_metadata_grch38.parquet"
    panel.to_parquet(out_meta, index=False)
    print(f"  sample_metadata_grch38.parquet: {len(panel)} samples")
    print(panel["super_pop"].value_counts().to_string())

    # PCA
    pcs = pd.read_csv(eigenvec, sep="\t")
    pcs = pcs.rename(columns={"#IID": "sample_id"})
    if "FID" in pcs.columns:
        pcs = pcs.drop(columns=["FID"])
    pcs = pcs.merge(panel, on="sample_id", how="left")
    out_pcs = root / "data/processed/ancestry_pcs_chr22.parquet"
    pcs.to_parquet(out_pcs, index=False)
    print(f"\n  ⚠️  CHR22 SMOKE TEST: ancestry_pcs_chr22.parquet ({len(pcs)} samples × {len(pcs.columns)} cols)")

    # Eigenvalues
    ev = pd.read_csv(eigenval, header=None, names=["eigenvalue"])
    ev["PC"] = [f"PC{i+1}" for i in range(len(ev))]
    ev["variance_explained"] = ev["eigenvalue"] / ev["eigenvalue"].sum()
    ev.to_parquet(processed / "chr22_pca_eigenvalues.parquet", index=False)
    print("\n  Top 5 PCs (chr22 only — inflated):")
    print(ev.head(5).to_string(index=False))


def main():
    print("=== FAIRGEN-Open: GRCh38 chr22 Smoke Test ===")
    print("⚠️  CHR22 ONLY — results not genome-wide\n")

    print("\n--- Step 1: VCF → pgen ---")
    step1_vcf_to_pgen()

    print("\n--- Step 2: QC PCA pgen (MAF>1%) ---")
    step2_qc_pca_pgen()

    print("\n--- Step 3: Scoring pgen (no MAF filter) ---")
    step3_score_pgen()

    print("\n--- Step 4: LD pruning ---")
    step4_ld_prune()

    print("\n--- Step 5: PCA (20 PCs) ---")
    step5_pca()

    print("\n--- Step 6: Save parquets ---")
    step6_save_parquets()

    print("\n=== chr22 smoke test complete ===")
    print("Next: run scripts/harmonize_grch38.py to harmonize 9 PGS → chr22 GRCh38")


if __name__ == "__main__":
    main()
