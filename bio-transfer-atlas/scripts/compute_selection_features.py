"""
Compute selection-turnover proxy features for a target chromosome.
"""
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data/features/selection"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def normalize_chr(value: str) -> str:
    s = str(value).strip()
    return s[3:] if s.lower().startswith("chr") else s


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22", help="Chromosome number/name (e.g. 22 or chr22)")
    return p.parse_args()


def maybe_write_legacy(chrom: str, src: Path, legacy_name: str) -> None:
    if chrom == "22":
        pd.read_parquet(src).to_parquet(OUT_DIR / legacy_name, index=False)


def resolve_af_path(chrom: str) -> Path:
    suffixed = ROOT / "data/features/af" / f"af_divergence_features.chr{chrom}.parquet"
    if suffixed.exists():
        return suffixed
    if chrom == "22":
        legacy = ROOT / "data/features/af/af_divergence_features.parquet"
        if legacy.exists():
            return legacy
    raise FileNotFoundError(
        f"Missing AF features for chr{chrom}. Run compute_af_features.py --chrom {chrom} first."
    )


def hudson_fst(p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
    num = (p1 - p2) ** 2
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return np.where(den > 0, np.clip(num / den, 0, 1), 0.0)


def main() -> None:
    args = parse_args()
    chrom = normalize_chr(args.chrom)
    af_path = resolve_af_path(chrom)
    af = pd.read_parquet(af_path)
    print(f"Loaded AF features for chr{chrom}: {len(af):,} variants")

    p = {sp: af[f"AF_{sp}"].clip(1e-6, 1 - 1e-6).values for sp in SUPERPOPS}
    sel = pd.DataFrame({"variant_id": af["variant_id"]})

    fst: dict[tuple[str, str], np.ndarray] = {}
    for a, b in combinations(SUPERPOPS, 2):
        f = hudson_fst(p[a], p[b])
        fst[(a, b)] = f
        fst[(b, a)] = f
        sel[f"FST_{a}_{b}"] = f

    def t_stat(a: str, b: str) -> np.ndarray:
        return -np.log(np.clip(1 - fst[(a, b)], 1e-6, 1.0))

    def pbs(focal: str, r1: str, r2: str) -> np.ndarray:
        return (t_stat(focal, r1) + t_stat(focal, r2) - t_stat(r1, r2)) / 2.0

    sel["PBS_EUR"] = pbs("EUR", "AFR", "EAS")
    sel["PBS_EAS"] = pbs("EAS", "EUR", "AFR")
    sel["PBS_AFR"] = pbs("AFR", "EUR", "EAS")
    sel["PBS_SAS"] = pbs("SAS", "EUR", "AFR")
    sel["PBS_AMR"] = pbs("AMR", "EUR", "AFR")

    pbs_cols = ["PBS_EUR", "PBS_EAS", "PBS_AFR", "PBS_SAS", "PBS_AMR"]
    sel["PBS_max"] = sel[pbs_cols].max(axis=1)
    sel["PBS_mean"] = sel[pbs_cols].mean(axis=1)
    sel["PBS_argmax_pop"] = sel[pbs_cols].idxmax(axis=1).str.replace("PBS_", "", regex=False)

    af_mat = np.vstack([p[sp] for sp in SUPERPOPS]).T
    sel["DAF_spread"] = af_mat.max(axis=1) - af_mat.min(axis=1)
    sel["DAF_std"] = af_mat.std(axis=1)

    out = OUT_DIR / f"selection_turnover_features.chr{chrom}.parquet"
    sel.to_parquet(out, index=False)
    maybe_write_legacy(chrom, out, "selection_turnover_features.parquet")
    print(f"Saved {out.name} ({len(sel):,} variants × {sel.shape[1] - 1} features)")

    nonzero = (sel["PBS_max"] > 0).mean() * 100
    print(f"\nGate chr{chrom}: {nonzero:.1f}% variants have a non-zero selection signal")


if __name__ == "__main__":
    main()
