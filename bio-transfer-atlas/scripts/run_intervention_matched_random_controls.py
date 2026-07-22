#!/usr/bin/env python3
"""
Matched random-n / random-mass Monte Carlo controls for interventions + LOSO.

Uses 1000G superpop AFs to compute expected ancestry mean scores under additive
weights (fast; enables 500+ draws). Observed MAD reduction still comes from
scored matrices when available.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from intervention_common import (  # noqa: E402
    GW_INT_ROOT,
    GW_SCORE_MATRIX,
    GW_TAG,
    INTERVENTION_MODES,
    PGS_IDS,
    PGS_TRAITS,
    SUPERPOPS,
    load_genomewide_weights,
)

TABLES = ROOT / "results/tables"
AF_PATH = ROOT / "data/features/af/1000g_af_by_superpop.parquet"
INT_SCORES = ROOT / "data/processed/scores_grch38_intervention_genomewide"
SEED = 719


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--n-draws", type=int, default=500)
    p.add_argument("--modes", default=",".join(m for m in INTERVENTION_MODES if not str(m).startswith("random")))
    p.add_argument("--pgs", default=",".join(PGS_IDS))
    p.add_argument("--intervention-root", default=str(GW_INT_ROOT))
    p.add_argument("--out", default=str(TABLES / "intervention_matched_random_controls.csv"))
    p.add_argument("--out-loso", default=str(TABLES / "intervention_loso_mad_by_mode.csv"))
    p.add_argument("--use-scored-mad", action="store_true", help="Also load scored matrices (slow)")
    return p.parse_args()


def mad_from_means(means: dict[str, float]) -> float:
    if "EUR" not in means or not np.isfinite(means["EUR"]):
        return float("nan")
    eur = means["EUR"]
    vals = [abs(means[sp] - eur) for sp in SUPERPOPS if sp != "EUR" and sp in means and np.isfinite(means[sp])]
    return float(np.mean(vals)) if vals else float("nan")


def contrib_matrix(weights: pd.DataFrame, af: pd.DataFrame | None = None) -> tuple[np.ndarray, list[str]]:
    """Per-variant expected score contribution by superpop: 2 * AF_effect * w."""
    if af is not None and "AF_EUR" not in weights.columns:
        m = weights.merge(af, on="variant_id", how="inner")
    else:
        m = weights
    sps = [sp for sp in SUPERPOPS if f"AF_{sp}" in m.columns]
    if m.empty or "effect_weight" not in m.columns or not sps:
        return np.zeros((0, 0)), sps
    ea = m["effect_allele"].astype(str).str.upper().to_numpy()
    ref = m["REF"].astype(str).str.upper().to_numpy() if "REF" in m.columns else np.array(["?"] * len(m))
    alt = m["ALT"].astype(str).str.upper().to_numpy() if "ALT" in m.columns else np.array(["?"] * len(m))
    w = pd.to_numeric(m["effect_weight"], errors="coerce").to_numpy(float)
    C = np.zeros((len(m), len(sps)), dtype=float)
    for j, sp in enumerate(sps):
        afv = pd.to_numeric(m[f"AF_{sp}"], errors="coerce").to_numpy(float)
        af_eff = np.where(ea == alt, afv, np.where(ea == ref, 1.0 - afv, np.nan))
        C[:, j] = 2.0 * af_eff * w
    return np.nan_to_num(C, nan=0.0), sps


def mad_from_contrib(C: np.ndarray, sps: list[str], keep: np.ndarray | None = None) -> float:
    if C.size == 0 or "EUR" not in sps:
        return float("nan")
    sub = C if keep is None else C[keep]
    means = {sp: float(sub[:, j].sum()) for j, sp in enumerate(sps)}
    return mad_from_means(means)


def expected_means(weights: pd.DataFrame, af: pd.DataFrame | None = None) -> dict[str, float]:
    """E[score] ≈ sum 2 * AF_effect * w after allele alignment."""
    C, sps = contrib_matrix(weights, af)
    if C.size == 0:
        return {sp: float("nan") for sp in SUPERPOPS}
    return {sp: float(C[:, j].sum()) for j, sp in enumerate(sps)}


def scored_mad(baseline: pd.DataFrame, edited: pd.DataFrame | None, pgs: str) -> tuple[float, float]:
    def _mad(df: pd.DataFrame) -> float:
        if pgs not in df.columns or "super_pop" not in df.columns:
            return float("nan")
        eur = df.loc[df["super_pop"] == "EUR", pgs].dropna()
        if eur.empty:
            return float("nan")
        m_eur = float(eur.mean())
        vals = []
        for sp, g in df.groupby("super_pop"):
            if sp == "EUR":
                continue
            vals.append(abs(float(g[pgs].mean()) - m_eur))
        return float(np.mean(vals)) if vals else float("nan")

    mad_b = _mad(baseline)
    mad_e = _mad(edited) if edited is not None else float("nan")
    return mad_b, mad_e


def match_mass_mask(abs_w: np.ndarray, target_removed: float, rng: np.random.Generator) -> np.ndarray:
    """Return boolean keep-mask after randomly removing ~target absolute weight mass."""
    n = len(abs_w)
    order = rng.permutation(n)
    cum = np.cumsum(abs_w[order])
    # drop until cumsum reaches target
    k = int(np.searchsorted(cum, target_removed, side="left")) + 1
    k = min(max(k, 1), n - 1)
    drop = np.zeros(n, dtype=bool)
    drop[order[:k]] = True
    return ~drop


def run_controls_for(
    pgs: str,
    mode: str,
    int_root: Path,
    o: pd.DataFrame,
    C: np.ndarray,
    sps: list[str],
    n_orig: int,
    w_orig: float,
    baseline_scores: pd.DataFrame | None,
    n_draws: int,
) -> dict:
    wpath = int_root / pgs / f"{mode}.tsv"
    if o.empty or not wpath.exists():
        return {"pgs_id": pgs, "mode": mode, "status": "missing_weights"}
    edited = pd.read_csv(wpath, sep="\t", usecols=lambda c: c in ("variant_id", "effect_weight"))
    if "effect_weight" not in edited.columns or "variant_id" not in edited.columns:
        return {"pgs_id": pgs, "mode": mode, "status": "bad_edited_weights"}

    edited = edited.dropna(subset=["effect_weight"]).copy()
    edit_ids = set(edited["variant_id"].astype(str))
    n_edit = len(edit_ids)
    n_removed = max(n_orig - n_edit, 0)
    w_edit = float(edited["effect_weight"].abs().sum())
    mass_removed = max(w_orig - w_edit, 0.0)
    variant_retention = n_edit / n_orig if n_orig else np.nan
    weight_mass_retention = w_edit / w_orig if w_orig else np.nan

    mad_b = mad_e = np.nan
    # Scored MAD is optional (slow parquet I/O); empirical_p uses AF-expected MAD.
    if baseline_scores is not None and getattr(run_controls_for, "_use_scored", False):
        epath = INT_SCORES / f"score_matrix_{mode}_{GW_TAG}.parquet"
        edited_sc = pd.read_parquet(epath) if epath.exists() else None
        mad_b, mad_e = scored_mad(baseline_scores, edited_sc, pgs)
    obs_reduction = mad_b - mad_e if np.isfinite(mad_b) and np.isfinite(mad_e) else np.nan

    n = len(o)
    if n < 50:
        return {
            "pgs_id": pgs,
            "mode": mode,
            "status": "insufficient_af_overlap",
            "n_af_overlap": n,
            "variant_retention": variant_retention,
            "weight_mass_retention": weight_mass_retention,
            "observed_mad_baseline": mad_b,
            "observed_mad_edited": mad_e,
            "observed_mad_reduction": obs_reduction,
        }

    base_exp_mad = mad_from_contrib(C, sps)
    edit_mask = o["variant_id"].astype(str).isin(edit_ids).to_numpy()
    edit_exp_mad = mad_from_contrib(C, sps, edit_mask)
    obs_exp_reduction = (
        base_exp_mad - edit_exp_mad if np.isfinite(base_exp_mad) and np.isfinite(edit_exp_mad) else np.nan
    )

    rng = np.random.default_rng(SEED + hash(pgs + mode) % 10_000)
    abs_w = o["effect_weight"].abs().to_numpy(float)
    n_rem = min(n_removed, n - 1) if n_removed > 0 else 0
    # Scale target mass to AF-overlap subset
    target_mass = mass_removed * (float(abs_w.sum()) / w_orig) if w_orig else 0.0

    red_n = np.empty(n_draws if n_rem > 0 else 0, dtype=float)
    red_m = np.empty(n_draws if target_mass > 0 else 0, dtype=float)
    for i in range(n_draws):
        if n_rem > 0:
            keep = np.ones(n, dtype=bool)
            keep[rng.choice(n, size=n_rem, replace=False)] = False
            red_n[i] = base_exp_mad - mad_from_contrib(C, sps, keep)
        if target_mass > 0:
            keep_m = match_mass_mask(abs_w, target_mass, rng)
            red_m[i] = base_exp_mad - mad_from_contrib(C, sps, keep_m)

    def emp_p(obs: float, draws: np.ndarray) -> float:
        if draws.size == 0 or not np.isfinite(obs):
            return float("nan")
        return float((1 + np.sum(draws >= obs)) / (1 + len(draws)))

    return {
        "pgs_id": pgs,
        "mode": mode,
        "trait": PGS_TRAITS.get(pgs, ""),
        "status": "ok",
        "n_variants_original": n_orig,
        "n_variants_edited": n_edit,
        "n_removed": n_removed,
        "variant_retention": variant_retention,
        "weight_mass_retention": weight_mass_retention,
        "mass_removed": mass_removed,
        "n_af_overlap": n,
        "n_draws": n_draws,
        "observed_mad_baseline": mad_b,
        "observed_mad_edited": mad_e,
        "observed_mad_reduction": obs_reduction,
        "expected_mad_baseline": base_exp_mad,
        "expected_mad_edited": edit_exp_mad,
        "expected_mad_reduction": obs_exp_reduction,
        "delta_vs_random_n": (
            obs_exp_reduction - float(red_n.mean()) if red_n.size and np.isfinite(obs_exp_reduction) else np.nan
        ),
        "delta_vs_random_mass": (
            obs_exp_reduction - float(red_m.mean()) if red_m.size and np.isfinite(obs_exp_reduction) else np.nan
        ),
        "empirical_p_n": emp_p(obs_exp_reduction, red_n),
        "empirical_p_mass": emp_p(obs_exp_reduction, red_m),
        "random_n_mad_reduction_mean": float(red_n.mean()) if red_n.size else np.nan,
        "random_mass_mad_reduction_mean": float(red_m.mean()) if red_m.size else np.nan,
        "claim": "ancestry_mean_separation_vs_matched_random_af_expected",
    }


def loso_table(control_df: pd.DataFrame) -> pd.DataFrame:
    mad = control_df[
        (control_df["status"] == "ok") & control_df["observed_mad_reduction"].notna()
    ].copy()
    if mad.empty:
        # fall back to expected
        mad = control_df[(control_df["status"] == "ok") & control_df["expected_mad_reduction"].notna()].copy()
        col = "expected_mad_reduction"
    else:
        col = "observed_mad_reduction"
    rows = []
    scores = sorted(mad["pgs_id"].unique())
    modes = sorted(mad["mode"].unique())
    for mode in modes:
        sub = mad[mad["mode"] == mode]
        rows.append(
            {
                "omit_pgs": "NONE",
                "mode": mode,
                "metric": col,
                "mean_mad_reduction": float(sub[col].mean()),
                "median_mad_reduction": float(sub[col].median()),
                "n_scores": len(sub),
            }
        )
        for omit in scores:
            rest = sub[sub["pgs_id"] != omit]
            if rest.empty:
                continue
            rows.append(
                {
                    "omit_pgs": omit,
                    "mode": mode,
                    "metric": col,
                    "mean_mad_reduction": float(rest[col].mean()),
                    "median_mad_reduction": float(rest[col].median()),
                    "n_scores": len(rest),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)
    af = pd.read_parquet(AF_PATH)
    af["variant_id"] = af["variant_id"].astype(str)
    af_cols = ["variant_id", "REF", "ALT"] + [f"AF_{sp}" for sp in SUPERPOPS if f"AF_{sp}" in af.columns]
    af = af[af_cols]
    baseline = None
    if args.use_scored_mad and Path(GW_SCORE_MATRIX).exists():
        baseline = pd.read_parquet(GW_SCORE_MATRIX)
        run_controls_for._use_scored = True  # type: ignore[attr-defined]
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    pgs_list = [p.strip() for p in args.pgs.split(",") if p.strip()]
    int_root = Path(args.intervention_root)
    rows = []
    for pgs in pgs_list:
        print(f"loading weights {pgs} ...", flush=True)
        try:
            orig = load_genomewide_weights(pgs)
        except FileNotFoundError as e:
            print(f"  skip {pgs}: {e}", flush=True)
            for mode in modes:
                rows.append({"pgs_id": pgs, "mode": mode, "status": "missing_weights"})
            continue
        if orig.empty:
            for mode in modes:
                rows.append({"pgs_id": pgs, "mode": mode, "status": "missing_weights"})
            continue
        orig = orig.dropna(subset=["effect_weight"]).copy()
        orig["variant_id"] = orig["variant_id"].astype(str)
        n_orig = int(orig["variant_id"].nunique())
        w_orig = float(orig["effect_weight"].abs().sum())
        o = orig.merge(af, on="variant_id", how="inner")
        C, sps = contrib_matrix(o)
        for mode in modes:
            print(f"controls {pgs} {mode} ...", flush=True)
            rows.append(
                run_controls_for(
                    pgs, mode, int_root, o, C, sps, n_orig, w_orig, baseline, args.n_draws
                )
            )
    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False, float_format="%.6g")
    loso = loso_table(out)
    loso.to_csv(args.out_loso, index=False, float_format="%.6g")
    print(f"Saved {args.out} ({len(out)} rows) and {args.out_loso}")


if __name__ == "__main__":
    main()
