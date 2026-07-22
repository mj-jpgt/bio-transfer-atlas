"""
Allele harmonization for PGS scoring files against a PLINK2 .pvar reference.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
AMBIGUOUS_PAIRS = {frozenset(["A", "T"]), frozenset(["C", "G"])}


def is_ambiguous(a1: str, a2: str) -> bool:
    return frozenset([a1.upper(), a2.upper()]) in AMBIGUOUS_PAIRS


def remove_ambiguous(df: pd.DataFrame) -> pd.DataFrame:
    mask = df.apply(
        lambda r: not is_ambiguous(r["effect_allele"], r["other_allele"]), axis=1
    )
    n_removed = (~mask).sum()
    if n_removed:
        logger.info(f"  Removed {n_removed} ambiguous A/T or C/G SNPs")
    return df[mask].copy()


def flip_alleles(row: dict[str, Any], ref_effect: str, ref_other: str) -> dict[str, Any]:
    r = dict(row)
    r["effect_allele"] = ref_effect
    r["other_allele"] = ref_other
    r["effect_weight"] = -float(row["effect_weight"])
    return r


def load_pgs_harmonized(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        lines = [l for l in f if not l.startswith("##")]
    from io import StringIO
    df = pd.read_csv(StringIO("".join(lines)), sep="\t", low_memory=False)
    df.columns = [c.lower() for c in df.columns]
    return df


def load_pvar(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path, sep="\t", comment="#",
        names=["chrom", "pos", "id", "ref", "alt"],
        usecols=["chrom", "pos", "id", "ref", "alt"],
    )
    df["pos"] = df["pos"].astype(int)
    return df


def harmonize(score_path: Path, pvar_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Match score variants to pvar, flip strands where needed, remove ambiguous.

    Returns
    -------
    harmonized : pd.DataFrame  — variants ready for PLINK scoring
    report     : pd.DataFrame  — match statistics
    """
    score = load_pgs_harmonized(score_path)
    pvar = load_pvar(pvar_path)

    score = remove_ambiguous(score)

    ref_lookup = {(str(r.chrom), int(r.pos)): r for r in pvar.itertuples()}

    matched, flipped, missing = [], 0, 0
    for _, row in score.iterrows():
        key = (str(row.get("chr_name", row.get("hm_chr", ""))),
               int(row.get("chr_position", row.get("hm_pos", 0))))
        ref = ref_lookup.get(key)
        if ref is None:
            missing += 1
            continue

        ea = str(row["effect_allele"]).upper()
        oa = str(row.get("other_allele", row.get("reference_allele", ""))).upper()
        ref_alleles = {ref.ref.upper(), ref.alt.upper()}

        if ea in ref_alleles and oa in ref_alleles:
            matched.append(row)
        elif (COMPLEMENT.get(ea, "") in ref_alleles and
              COMPLEMENT.get(oa, "") in ref_alleles):
            flipped += 1
            row = row.copy()
            row["effect_allele"] = COMPLEMENT[ea]
            row["other_allele"] = COMPLEMENT[oa]
            matched.append(row)
        else:
            missing += 1

    harmonized = pd.DataFrame(matched) if matched else pd.DataFrame(columns=score.columns)

    report = pd.DataFrame([{
        "n_variants_original": len(score),
        "n_variants_matched": len(harmonized),
        "n_flipped": flipped,
        "n_missing": missing,
        "match_rate": len(harmonized) / max(len(score), 1),
    }])

    return harmonized, report
