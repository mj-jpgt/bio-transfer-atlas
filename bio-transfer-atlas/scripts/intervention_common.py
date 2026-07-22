"""Shared constants and helpers for Phase 18 intervention pilot."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from gwas_utils import normalize_chr

ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(ROOT / "scripts"))
from bta_plink import plink2_bin  # noqa: E402

PLINK2 = plink2_bin()
SEED = 719

PGS_TRAITS = {
    "PGS000018": "T2D",
    "PGS004696": "CAD",
    "PGS004698": "CAD",
    "PGS003897": "BMI",
    "PGS002853": "LDL",
    "PGS002858": "LDL",
    "PGS003092": "BMI",
    "PGS000014": "CAD",
    "PGS004840": "T2D",
    # M4 expansions
    "PGS000191": "WBC",  # Duffy positive control
    "PGS004133": "RA",   # MHC stress-test
    "PGS001288": "IBD",  # MHC stress-test
}
INTERVENTION_MODES_EXTRA = ["duffy_gate"]  # ACKR1/rs2814778 window drop
PGS_IDS = list(PGS_TRAITS.keys())
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

INTERVENTION_MODES = [
    "filter_5",
    "filter_10",
    "filter_20",
    "reweight_linear",
    "reweight_exp",
    "flag",
    "random",
    "fst",
    "ld",
    "maf",
    "duffy_gate",
]
# Duffy-null rs2814778 (GRCh38)
# Classic Duffy-null allele is C (FY*01N.01 / rs2814778-C), common in AFR, rare in EUR.
DUFFY_CHR = "1"
DUFFY_POS = 159_204_893
DUFFY_RSID = "rs2814778"
DUFFY_NULL_ALLELE = "C"  # allele counted as dosage for "null" biology
DUFFY_WINDOW_BP = 200_000
FILTER_FRACS = {"filter_5": 0.05, "filter_10": 0.10, "filter_20": 0.20}
NEG_CONTROL_FRAC = 0.10
REWEIGHT_EXP_K = 3.0
GW_CHROMS = [str(i) for i in range(1, 23)]
INTERIM = ROOT / "data/interim/1000g_grch38"
MDIR = ROOT / "data/modeling"
GW_TAG = "genomewide"
GW_MASTER = MDIR / f"master_variant_table_genomewide_{GW_TAG}.parquet"
GW_MODEL = MDIR / f"portability_model.{GW_TAG}.joblib"
GW_PREDS = MDIR / f"variant_portability_predictions.{GW_TAG}.parquet"
GW_INT_ROOT = ROOT / "data/processed/pgs_grch38_intervention_genomewide"
GW_SCORE_MATRIX = ROOT / "data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet"


def is_chr22_variant(variant_id: str) -> bool:
    return normalize_chr(str(variant_id).split(":")[0]) == "22"


def filter_chr22(df: pd.DataFrame, col: str = "variant_id") -> pd.DataFrame:
    return df[df[col].map(is_chr22_variant)].copy()


def harmonized_path(pgs_id: str, pgs_root: Path | None = None) -> Path:
    root = pgs_root or (ROOT / "data/processed/pgs_grch38")
    return root / pgs_id / f"{pgs_id}.harmonized.tsv"


def load_chr22_weights(pgs_id: str, pgs_root: Path | None = None) -> pd.DataFrame:
    return load_pgs_weights(pgs_id, chroms=["22"], pgs_root=pgs_root)


def variant_chrom(variant_id: str) -> str:
    return normalize_chr(str(variant_id).split(":")[0])


def pfile_for_chrom(chrom: str) -> Path:
    return INTERIM / f"chr{chrom}.score"


def available_score_chroms(chroms: list[str] | None = None) -> list[str]:
    wanted = chroms or GW_CHROMS
    out = []
    for c in wanted:
        if Path(str(pfile_for_chrom(c)) + ".pgen").exists():
            out.append(str(int(c)))
    return out


def load_pgs_weights(
    pgs_id: str,
    chroms: list[str] | None = None,
    pgs_root: Path | None = None,
) -> pd.DataFrame:
    """Load harmonized PGS weights, optionally filtered to chromosome list."""
    path = harmonized_path(pgs_id, pgs_root)
    chrom_set = {normalize_chr(c) for c in chroms} if chroms else None
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        chunksize=250_000,
        dtype={"variant_id": str, "effect_allele": str},
    ):
        chunk["effect_weight"] = pd.to_numeric(chunk["effect_weight"], errors="coerce")
        if chrom_set is not None:
            sub = chunk[chunk["variant_id"].map(lambda v: variant_chrom(v) in chrom_set)]
        else:
            sub = chunk
        if not sub.empty:
            chunks.append(sub)
    if not chunks:
        return pd.DataFrame(columns=["variant_id", "effect_allele", "effect_weight"])
    return pd.concat(chunks, ignore_index=True)


def load_genomewide_weights(pgs_id: str, pgs_root: Path | None = None) -> pd.DataFrame:
    return load_pgs_weights(pgs_id, chroms=None, pgs_root=pgs_root)
