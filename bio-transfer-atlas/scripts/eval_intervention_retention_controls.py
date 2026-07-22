#!/usr/bin/env python3
"""
Intervention retention / variance / correlation metrics + retention–MAD curve.

Reframes outcomes as ancestry mean-separation (not portability/accuracy).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results/tables"
BASE = ROOT / "data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet"
INT = ROOT / "data/processed/scores_grch38_intervention_genomewide"
INT_ROOT = ROOT / "data/processed/pgs_grch38_intervention_genomewide"


def mad(df: pd.DataFrame, col: str) -> float:
    eur = df.loc[df["super_pop"] == "EUR", col]
    if eur.empty:
        return float("nan")
    m = float(eur.mean())
    vals = []
    for sp, g in df.groupby("super_pop"):
        if sp == "EUR":
            continue
        vals.append(abs(float(g[col].mean()) - m))
    return float(np.mean(vals)) if vals else float("nan")


def main() -> None:
    if not BASE.exists():
        print("missing baseline scores")
        return
    base = pd.read_parquet(BASE).dropna(subset=["super_pop"])
    pgs_cols = [c for c in base.columns if str(c).startswith("PGS")]
    # Keep runtime bounded: prioritize known metabolic/WBC/autoimmune scores
    prefer = [
        "PGS000018",
        "PGS000011",
        "PGS000027",
        "PGS000191",
        "PGS004133",
        "PGS001288",
    ]
    ordered = [c for c in prefer if c in pgs_cols] + [c for c in pgs_cols if c not in prefer]
    pgs_cols = ordered[:8]
    rows = []
    curve = []
    modes = []
    if INT.exists():
        modes = sorted(
            p.name.replace("score_matrix_", "").replace("_genomewide.parquet", "")
            for p in INT.glob("score_matrix_*_genomewide.parquet")
        )
    # Prefer informative modes + random controls; cap for runtime
    prefer_modes = [
        "flag",
        "fst",
        "maf",
        "filter_10",
        "random",
        "random_10",
        "duffy_gate",
        "reweight_linear",
    ]
    modes = [m for m in prefer_modes if m in modes] + [
        m for m in modes if m not in prefer_modes and m.startswith("random")
    ]
    modes = modes[:10]
    for mode in modes:
        path = INT / f"score_matrix_{mode}_genomewide.parquet"
        if not path.exists():
            continue
        sc = pd.read_parquet(path).dropna(subset=["super_pop"])
        for pgs in pgs_cols:
            if pgs not in sc.columns or pgs not in base.columns:
                continue
            # Align samples
            m = base[["sample_id", "super_pop", pgs]].merge(
                sc[["sample_id", pgs]], on="sample_id", suffixes=("_base", "_edit")
            )
            m = m.rename(columns={f"{pgs}_base": "base", f"{pgs}_edit": "edit"})
            mad_b = mad(m.rename(columns={"base": pgs}), pgs)
            mad_e = mad(m.rename(columns={"edit": pgs}), pgs)
            # Variance retained
            var_rows = {}
            corr_rows = {}
            for sp, g in m.groupby("super_pop"):
                vb = float(g["base"].var())
                ve = float(g["edit"].var())
                var_rows[sp] = ve / vb if vb > 0 else np.nan
                if len(g) > 5:
                    corr_rows[sp] = float(np.corrcoef(g["base"], g["edit"])[0, 1])
            # Weight mass from TSVs if present
            w_ret = np.nan
            n_ret = np.nan
            wpath = INT_ROOT / pgs / f"{mode}.tsv"
            wbase = None
            # try harmonized as baseline weight mass
            harm = ROOT / "data/processed/pgs_grch38" / pgs / f"{pgs}.harmonized.tsv"
            if harm.exists():
                wb = pd.read_csv(harm, sep="\t")
                wsum = wb["effect_weight"].abs().sum() if "effect_weight" in wb.columns else np.nan
                n0 = len(wb)
                if wpath.exists():
                    we = pd.read_csv(wpath, sep="\t")
                    w_ret = float(we["effect_weight"].abs().sum() / wsum) if wsum else np.nan
                    n_ret = float(len(we) / n0) if n0 else np.nan
                elif mode == "flag":
                    w_ret, n_ret = 1.0, 1.0
            rows.append(
                {
                    "pgs_id": pgs,
                    "mode": mode,
                    "metric": "mean_abs_delta_EUR",
                    "value": mad_e,
                    "baseline_value": mad_b,
                    "reduction": mad_b - mad_e if np.isfinite(mad_b) and np.isfinite(mad_e) else np.nan,
                    "claim": "ancestry_mean_separation_not_portability",
                }
            )
            rows.append(
                {
                    "pgs_id": pgs,
                    "mode": mode,
                    "metric": "weight_mass_retained",
                    "value": w_ret,
                }
            )
            rows.append(
                {
                    "pgs_id": pgs,
                    "mode": mode,
                    "metric": "variant_frac_retained",
                    "value": n_ret,
                }
            )
            for sp, v in var_rows.items():
                rows.append(
                    {
                        "pgs_id": pgs,
                        "mode": mode,
                        "metric": f"score_var_retained_{sp}",
                        "value": v,
                    }
                )
            for sp, v in corr_rows.items():
                rows.append(
                    {
                        "pgs_id": pgs,
                        "mode": mode,
                        "metric": f"corr_edited_original_{sp}",
                        "value": v,
                    }
                )
            curve.append(
                {
                    "pgs_id": pgs,
                    "mode": mode,
                    "weight_mass_retained": w_ret,
                    "variant_frac_retained": n_ret,
                    "mad_reduction": mad_b - mad_e if np.isfinite(mad_b) and np.isfinite(mad_e) else np.nan,
                    "mad": mad_e,
                }
            )
    out = TABLES / "intervention_retention_variance_metrics.csv"
    TABLES.mkdir(parents=True, exist_ok=True)
    metrics = pd.DataFrame(rows)
    # Matched random-n / random-mass controls: compare structured modes to random_* modes
    # with closest variant_frac or weight_mass retained.
    control_rows = []
    if not metrics.empty and "metric" in metrics.columns:
        mad_df = metrics[metrics["metric"] == "mean_abs_delta_EUR"].copy()
        nret = metrics[metrics["metric"] == "variant_frac_retained"][
            ["pgs_id", "mode", "value"]
        ].rename(columns={"value": "n_ret"})
        wret = metrics[metrics["metric"] == "weight_mass_retained"][
            ["pgs_id", "mode", "value"]
        ].rename(columns={"value": "w_ret"})
        mad_df = mad_df.merge(nret, on=["pgs_id", "mode"], how="left").merge(
            wret, on=["pgs_id", "mode"], how="left"
        )
        rand_modes = [m for m in mad_df["mode"].unique() if str(m).startswith("random")]
        structured = [m for m in mad_df["mode"].unique() if m not in rand_modes and m != "flag"]
        for pgs, g in mad_df.groupby("pgs_id"):
            rsub = g[g["mode"].isin(rand_modes)]
            for mode in structured:
                row = g[g["mode"] == mode]
                if row.empty:
                    continue
                row = row.iloc[0]
                best_n, best_w = None, None
                if len(rsub) and pd.notna(row.get("n_ret")):
                    best_n = rsub.iloc[(rsub["n_ret"] - row["n_ret"]).abs().argmin()]
                if len(rsub) and pd.notna(row.get("w_ret")):
                    best_w = rsub.iloc[(rsub["w_ret"] - row["w_ret"]).abs().argmin()]
                control_rows.append(
                    {
                        "pgs_id": pgs,
                        "mode": mode,
                        "metric": "mad_vs_matched_random_n",
                        "value": float(row["value"])
                        - (float(best_n["value"]) if best_n is not None else float("nan")),
                        "matched_random_mode": best_n["mode"] if best_n is not None else None,
                        "claim": "ancestry_mean_separation_vs_matched_random",
                    }
                )
                control_rows.append(
                    {
                        "pgs_id": pgs,
                        "mode": mode,
                        "metric": "mad_vs_matched_random_mass",
                        "value": float(row["value"])
                        - (float(best_w["value"]) if best_w is not None else float("nan")),
                        "matched_random_mode": best_w["mode"] if best_w is not None else None,
                        "claim": "ancestry_mean_separation_vs_matched_random",
                    }
                )
    if control_rows:
        metrics = pd.concat([metrics, pd.DataFrame(control_rows)], ignore_index=True)
    metrics.to_csv(out, index=False, float_format="%.6g")
    curvep = TABLES / "intervention_retention_mad_curve.csv"
    pd.DataFrame(curve).to_csv(curvep, index=False, float_format="%.6g")
    # Figure
    try:
        import matplotlib.pyplot as plt

        c = pd.DataFrame(curve).dropna(subset=["weight_mass_retained", "mad_reduction"])
        fig, ax = plt.subplots(figsize=(6.5, 4))
        for mode, g in c.groupby("mode"):
            ax.scatter(g["weight_mass_retained"], g["mad_reduction"], label=mode, alpha=0.8)
        ax.axhline(0, color="gray", ls="--", lw=0.8)
        ax.set_xlabel("Weight mass retained")
        ax.set_ylabel("MAD reduction (mean separation)")
        ax.set_title("Ancestry mean-separation vs retained information")
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(ROOT / "results/figures/fig_intervention_retention_curve.png", dpi=140)
        plt.close(fig)
    except Exception as e:
        print("figure skip", e)
    print(f"Saved {out} and {curvep}")


if __name__ == "__main__":
    main()
