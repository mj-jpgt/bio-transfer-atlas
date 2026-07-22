#!/usr/bin/env python3
"""
Paper figures for Biological Transferability Atlas (M5).

Writes PNGs under results/figures/ from existing results/tables CSVs.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results/tables"
FIGS = ROOT / "results/figures"
FIGS.mkdir(parents=True, exist_ok=True)


def _plt():
    import matplotlib.pyplot as plt

    plt.rcParams.update({"figure.dpi": 140, "font.size": 10})
    return plt


def fig_roc_ldblock() -> None:
    plt = _plt()
    path = TABLES / "ablation_ldblock_and_baselines_genomewide.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    sub = df[df["split"] == "split_ld_block"].copy()
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(6.5, 4))
    # AUROC bars (full ROC curves need predictions; report AUROC bars)
    order = sub.sort_values("AUROC", ascending=False)
    ax.barh(order["feature_group"], order["AUROC"], color="#2c5f7c")
    ax.axvline(0.5, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("AUROC (LD-block CV)")
    ax.set_title("Portability failure prediction by feature group")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_roc_ldblock.png")
    plt.close(fig)


def fig_feature_importance() -> None:
    plt = _plt()
    path = TABLES / "shap_mechanism_attribution_genomewide.csv"
    if not path.exists():
        # try gain table
        cands = list(TABLES.glob("*shap*.csv")) + list(TABLES.glob("*feature_importance*.csv"))
        if not cands:
            return
        path = cands[0]
    df = pd.read_csv(path)
    # flexible columns
    feat_col = next((c for c in ["feature", "feat", "name"] if c in df.columns), None)
    imp_col = next(
        (c for c in ["mean_abs_shap", "importance", "gain", "mean(|SHAP|)", "shap"] if c in df.columns),
        None,
    )
    if feat_col is None or imp_col is None:
        return
    top = df.nlargest(15, imp_col)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.barh(top[feat_col], top[imp_col], color="#3d6b51")
    ax.set_xlabel(imp_col)
    ax.set_title("Feature attribution (AF-dominant)")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_feature_importance.png")
    plt.close(fig)


def fig_duffy_shap() -> None:
    plt = _plt()
    strata = TABLES / "duffy_wbc_genotype_strata.csv"
    if not strata.exists():
        return
    df = pd.read_csv(strata)
    dose_col = "duffy_null_gt" if "duffy_null_gt" in df.columns else "duffy_gt"
    # Numeric genotype doses only (skip dominant labels / underpowered optional)
    num = df.copy()
    num["_dose"] = pd.to_numeric(num[dose_col], errors="coerce")
    num = num[num["_dose"].notna()]
    if "underpowered" in num.columns:
        num = num[~num["underpowered"].astype(bool)]
    fig, ax = plt.subplots(figsize=(6, 4))
    for sp, color in [("AFR", "#c45c26"), ("EUR", "#2c5f7c")]:
        sub = num[num["super_pop"] == sp]
        if sub.empty:
            continue
        yerr = None
        if {"mean_lo", "mean_hi"}.issubset(sub.columns):
            yerr = [
                sub["mean_wbc_pgs"] - sub["mean_lo"],
                sub["mean_hi"] - sub["mean_wbc_pgs"],
            ]
        elif "sd_wbc_pgs" in sub.columns:
            yerr = sub["sd_wbc_pgs"] / np.sqrt(sub["n"].clip(lower=1))
        ax.errorbar(
            sub["_dose"],
            sub["mean_wbc_pgs"],
            yerr=yerr,
            marker="o",
            label=sp,
            color=color,
            capsize=3,
        )
    ax.set_xlabel("Copies of Duffy-null allele (C)")
    ax.set_ylabel("Mean WBC PGS")
    ax.set_title("Duffy positive control (allele-audited)")
    ax.legend()
    ax.set_xticks([0, 1, 2])
    fig.tight_layout()
    fig.savefig(FIGS / "fig_duffy_shap.png")
    plt.close(fig)


def fig_intervention_mad() -> None:
    plt = _plt()
    path = TABLES / "intervention_results.genomewide.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    mad = df[df["metric"].astype(str).str.contains("mean_abs|MAD|mad", case=False, na=False)].copy()
    if mad.empty and "value" in df.columns:
        mad = df.copy()
    # Prefer mean_abs_delta_EUR
    if "metric" in mad.columns:
        prefer = mad[mad["metric"] == "mean_abs_delta_EUR"]
        if not prefer.empty:
            mad = prefer
    # reduction vs baseline if columns exist
    if "mode" not in mad.columns:
        return
    # Aggregate mean across PGS
    val_col = next((c for c in ["value", "mad", "mean_abs_delta", "estimate"] if c in mad.columns), None)
    if val_col is None:
        return
    g = mad.groupby("mode", as_index=False)[val_col].mean()
    # Also CI if present
    lo_col = next((c for c in ["ci_lo", "lo", "value_lo"] if c in mad.columns), None)
    hi_col = next((c for c in ["ci_hi", "hi", "value_hi"] if c in mad.columns), None)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    order = g.sort_values(val_col)
    yerr = None
    if lo_col and hi_col and lo_col in mad.columns:
        # approximate: mean of bounds
        lo = mad.groupby("mode")[lo_col].mean().reindex(order["mode"])
        hi = mad.groupby("mode")[hi_col].mean().reindex(order["mode"])
        yerr = np.vstack([order[val_col].to_numpy() - lo.to_numpy(), hi.to_numpy() - order[val_col].to_numpy()])
    ax.barh(order["mode"], order[val_col], xerr=yerr, color="#7a4e2d", capsize=2)
    ax.set_xlabel("Mean |Δ EUR| (MAD proxy)")
    ax.set_title("Intervention bake-off (genome-wide)")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_intervention_mad.png")
    plt.close(fig)


def fig_ancestry_means_fst() -> None:
    plt = _plt()
    # Use intervention summary or score shifts if available
    path = TABLES / "intervention_summary.genomewide.txt"
    # Prefer parquet score matrix means if present — else skip gracefully
    shifts = TABLES / "score_shifts_genomewide.csv"
    if not shifts.exists():
        cands = list(TABLES.glob("*score_shift*")) + list(TABLES.glob("*ancestry_mean*"))
        if not cands:
            return
        shifts = cands[0]
    if shifts.suffix == ".txt":
        return
    df = pd.read_csv(shifts) if shifts.suffix == ".csv" else pd.read_parquet(shifts)
    # flexible plot: if has ancestry + mean columns
    anc = next((c for c in ["super_pop", "ancestry", "population"] if c in df.columns), None)
    mean_c = next((c for c in ["mean_score", "mean", "value"] if c in df.columns), None)
    if anc is None or mean_c is None:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    for mode, sub in df.groupby(df.get("mode", pd.Series(["baseline"] * len(df)))):
        ax.plot(sub[anc], sub[mean_c], marker="o", label=str(mode))
    ax.set_title("Ancestry means (pre/post intervention)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_ancestry_means_fst.png")
    plt.close(fig)


def fig_finemap_tiers() -> None:
    plt = _plt()
    path = TABLES / "ablation_finemap_tiers.csv"
    if not path.exists():
        path = TABLES / "ablation_finemap_tiers_zlead.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    if "tier" in df.columns and "AUROC" in df.columns:
        for gname, sub in df.groupby(df.get("feature_group", pd.Series(["AF_LD_SEL"] * len(df)))):
            ax.bar(
                [f"{t}\n{gname}" for t in sub["tier"]],
                sub["AUROC"],
                alpha=0.8,
                label=str(gname),
            )
    ax.set_ylabel("AUROC")
    ax.set_title("Fine-mapping tier portability")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_finemap_tiers.png")
    plt.close(fig)


def fig_rg_ablation() -> None:
    """Peer-contest AUROC bars (excludes trait-constant Z diagnostic)."""
    plt = _plt()
    path = TABLES / "ablation_ldblock_peer_contest.csv"
    if not path.exists():
        path = TABLES / "ablation_ldblock_and_baselines_genomewide.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    sub = df[df["split"] == "split_ld_block"] if "split" in df.columns else df
    keep = sub[
        sub["feature_group"].isin(["AF_LD_SEL", "FST", "POP_DISTANCE", "RG_PROXY"])
        & ~sub["feature_group"].astype(str).str.contains("TRAIT_CONSTANT|RG_REAL", na=False)
    ]
    if keep.empty:
        return
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    colors = ["#2c5f7c", "#888", "#555", "#7a4e2d"][: len(keep)]
    ax.bar(keep["feature_group"], keep["AUROC"], color=colors)
    ax.set_ylabel("AUROC")
    ax.set_title("Variant-scale peer contest (LD-block CV)")
    ax.set_ylim(0.5, max(0.75, keep["AUROC"].max() + 0.05))
    fig.tight_layout()
    fig.savefig(FIGS / "fig_rg_ablation.png")
    plt.close(fig)


def main() -> None:
    fig_roc_ldblock()
    fig_feature_importance()
    fig_duffy_shap()
    fig_intervention_mad()
    fig_ancestry_means_fst()
    fig_finemap_tiers()
    fig_rg_ablation()
    print("Wrote figures under", FIGS)
    for p in sorted(FIGS.glob("fig_*.png")):
        print(" ", p.name)


if __name__ == "__main__":
    main()
