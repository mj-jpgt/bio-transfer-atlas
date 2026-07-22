"""
FAIRGEN-Open Stage 1-3 Data Downloads
======================================
Downloads:
  - 1000G GRCh38 high-coverage chr22 VCF + index (smoke test)
  - 1000G sample metadata panel (3202 samples)
  - PGS Catalog bulk metadata CSVs
  - 9 selected PGS GRCh38 harmonized scoring files
  - Reactome pathway files
  - LiftOver chain files (hg19<->hg38)

Run genome-wide download separately after smoke test passes.
Usage: python scripts/download_stage1_data.py [--genome-wide]
"""
import argparse
import time
from pathlib import Path

import requests
from tqdm import tqdm

root = Path(__file__).resolve().parents[1]

# ── helpers ────────────────────────────────────────────────────────────────
def download(url: str, dest: Path, desc: str = None) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  SKIP (exists): {dest.name}")
        return True
    tmp = dest.with_suffix(dest.suffix + ".partial")
    resume = tmp.stat().st_size if tmp.exists() else 0
    headers = {"Range": f"bytes={resume}-"} if resume else {}
    try:
        r = requests.get(url, headers=headers, stream=True, timeout=60)
        if r.status_code == 416:
            tmp.rename(dest); return True
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0)) + resume
        mode = "ab" if resume else "wb"
        with open(tmp, mode) as f, tqdm(
            desc=desc or dest.name, total=total, initial=resume,
            unit="B", unit_scale=True, unit_divisor=1024, leave=False
        ) as bar:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk); bar.update(len(chunk))
        tmp.rename(dest)
        print(f"  OK  {dest.name}  ({dest.stat().st_size/1e6:.1f} MB)")
        return True
    except Exception as e:
        print(f"  FAIL {dest.name}: {e}")
        return False

# ── 1. 1000G GRCh38 high-coverage VCFs ─────────────────────────────────────
HC_BASE = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
    "1kGP_high_coverage_Illumina.chr{c}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)
VCF_DIR = root / "data/raw/1000g/vcf_grch38"

def download_vcfs(chrs):
    print(f"\n=== 1000G GRCh38 VCFs: chr{chrs[0]}–chr{chrs[-1]} ===")
    for c in chrs:
        vcf_url = HC_BASE.format(c=c)
        tbi_url = vcf_url + ".tbi"
        download(vcf_url, VCF_DIR / f"chr{c}.vcf.gz", f"chr{c} VCF")
        download(tbi_url, VCF_DIR / f"chr{c}.vcf.gz.tbi", f"chr{c} .tbi")
        time.sleep(0.3)

# ── 2. Sample metadata ──────────────────────────────────────────────────────
META_DIR = root / "data/raw/1000g/metadata"

def download_metadata():
    print("\n=== 1000G Sample Metadata ===")
    urls = [
        ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
         "integrated_call_samples_v3.20130502.ALL.panel",
         "integrated_call_samples_v3.20130502.ALL.panel"),
        ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
         "1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt",
         "1000g_3202_samples_ped_population.txt"),
        ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
         "1000G_2504_high_coverage/1000G_2504_and_698_related_high_coverage.sequence.index",
         "1000G_2504_and_698_related_high_coverage.sequence.index"),
    ]
    for url, fname in urls:
        download(url, META_DIR / fname)

# ── 3. PGS Catalog bulk metadata ────────────────────────────────────────────
PGS_META_DIR = root / "data/raw/pgs_catalog/metadata"
PGS_META_URLS = {
    "pgs_all_metadata_scores.csv":
        "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/pgs_all_metadata_scores.csv",
    "pgs_all_metadata_publications.csv":
        "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/pgs_all_metadata_publications.csv",
    "pgs_all_metadata_evaluation_sample_sets.csv":
        "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/pgs_all_metadata_evaluation_sample_sets.csv",
    "pgs_all_metadata_score_development_samples.csv":
        "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/pgs_all_metadata_score_development_samples.csv",
}

def download_pgs_metadata():
    print("\n=== PGS Catalog Bulk Metadata ===")
    for fname, url in PGS_META_URLS.items():
        download(url, PGS_META_DIR / fname)

# ── 4. Selected PGS GRCh38 scoring files ────────────────────────────────────
PGS_SCORES_DIR = root / "data/raw/pgs_catalog/scores"
SELECTED_PGS = [
    "PGS000018",  # T2D - EUR
    "PGS004696",  # CAD - multi-ancestry
    "PGS004698",  # CAD - multi-ancestry
    "PGS003897",  # BMI - multi-ancestry
    "PGS002853",  # LDL - EUR
    "PGS002858",  # LDL - EAS
    "PGS003092",  # BMI - EAS
    "PGS000014",  # CAD - EUR
    "PGS004840",  # T2D - SAS
    "PGS000191",  # WBC - Duffy control
    "PGS004133",  # RA - MHC stress-test
    "PGS001288",  # IBD - MHC stress-test
]
PGS_FTP = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs}/ScoringFiles/Harmonized/{pgs}_hmPOS_GRCh38.txt.gz"

def download_pgs_scores():
    print("\n=== PGS GRCh38 Harmonized Scoring Files ===")
    for pgs in SELECTED_PGS:
        url = PGS_FTP.format(pgs=pgs)
        dest = PGS_SCORES_DIR / pgs / f"{pgs}_hmPOS_GRCh38.txt.gz"
        download(url, dest, f"{pgs} GRCh38")

# ── 5. Reactome ─────────────────────────────────────────────────────────────
REACTOME_DIR = root / "data/raw/reactome"
REACTOME_URLS = {
    "ReactomePathways.txt":
        "https://reactome.org/download/current/ReactomePathways.txt",
    "ReactomePathwaysRelation.txt":
        "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
    "NCBI2Reactome.txt":
        "https://reactome.org/download/current/NCBI2Reactome.txt",
    "Ensembl2Reactome.txt":
        "https://reactome.org/download/current/Ensembl2Reactome.txt",
    "UniProt2Reactome.txt":
        "https://reactome.org/download/current/UniProt2Reactome.txt",
}

def download_reactome():
    print("\n=== Reactome Pathway Files ===")
    for fname, url in REACTOME_URLS.items():
        download(url, REACTOME_DIR / fname)

# ── 6. LiftOver chain files ──────────────────────────────────────────────────
REF_DIR = root / "data/raw/reference"
CHAIN_URLS = {
    "hg19ToHg38.over.chain.gz":
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    "hg38ToHg19.over.chain.gz":
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
}

def download_reference():
    print("\n=== LiftOver Chain Files ===")
    for fname, url in CHAIN_URLS.items():
        download(url, REF_DIR / fname)

# ── main ────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome-wide", action="store_true",
                        help="Download all chr1-22 VCFs (26+ GB). Default: chr22 only.")
    args = parser.parse_args()

    chrs = list(range(1, 23)) if args.genome_wide else [22]
    if not args.genome_wide:
        print("*** CHR22 SMOKE TEST ONLY — pass --genome-wide to download all autosomes ***")

    download_vcfs(chrs)
    download_metadata()
    download_pgs_metadata()
    download_pgs_scores()
    download_reactome()
    download_reference()

    print("\n=== Download complete ===")

if __name__ == "__main__":
    main()
