"""
FAIRGEN-Open Stage 7 (cont.): Atlas figures
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.use("Agg")

root = Path(__file__).resolve().parents[1]
tables = root / "results/tables"
figs = root / "results/figures"
figs.mkdir(parents=True, exist_ok=True)

SUPERPOPS = ["AFR", "AMR", "SAS", "EAS", "EUR"]
TRAIT_COLORS = {"CAD": "#d62728", "T2D": "#1f77b4", "BMI": "#2ca02c", "LDL": "#ff7f0e"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Render atlas figures from score-shift tables.")
    p.add_argument("--tag", default="genomewide", help="Table suffix, e.g. genomewide or empty for chr22")
    return p.parse_args()


def resolve_table(name: str, tag: str) -> Path:
    if tag:
        tagged = tables / f"{name}_{tag}.parquet"
        if tagged.exists():
            return tagged
    legacy = tables / f"{name}.parquet"
    if legacy.exists():
        return legacy
    raise FileNotFoundError(f"Missing {name} table (tag={tag!r})")


def main() -> None:
    args = parse_args()
    tag = args.tag
    scope = "genome-wide" if tag == "genomewide" else (tag or "chr22")

    sp = pd.read_parquet(resolve_table("score_shifts_superpop", tag))
    ds = pd.read_parquet(resolve_table("distance_sensitivity", tag))
    sp["label"] = sp["pgs_id"] + " (" + sp["trait"] + ")"

    pivot = sp.pivot(index="label", columns="super_pop", values="delta_EUR")
    pivot = pivot.reindex(columns=SUPERPOPS).sort_index()

    fig, ax = plt.subplots(figsize=(7, 6))
    vmax = np.nanmax(np.abs(pivot.values))
    sns.heatmap(
        pivot, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
        annot=True, fmt=".2f", linewidths=0.5, ax=ax,
        cbar_kws={"label": "(mean_pop - mean_EUR) / SD_EUR"},
    )
    ax.set_title(f"PRS Score-Shift Atlas ({scope})\nEUR-centred standardized shift", fontsize=12)
    ax.set_xlabel("Super-population")
    plt.tight_layout()
    suffix = f"_{tag}" if tag else ""
    fig.savefig(figs / f"fig_atlas_heatmap{suffix}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved fig_atlas_heatmap{suffix}.png")

    ds_sorted = ds.sort_values("slope", ascending=True)
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = [TRAIT_COLORS.get(t, "#888") for t in ds_sorted["trait"]]
    y = np.arange(len(ds_sorted))
    ax.barh(
        y, ds_sorted["slope"], color=colors,
        xerr=[ds_sorted["slope"] - ds_sorted["slope_lo"], ds_sorted["slope_hi"] - ds_sorted["slope"]],
        capsize=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(ds_sorted["pgs_id"] + " (" + ds_sorted["trait"] + ")")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Distance-sensitivity slope  (|shift| per unit PC distance)")
    ax.set_title(f"Genetic-distance sensitivity by PGS ({scope})\n95% bootstrap CI", fontsize=12)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in TRAIT_COLORS.values()]
    ax.legend(handles, TRAIT_COLORS.keys(), title="Trait", loc="lower right")
    plt.tight_layout()
    fig.savefig(figs / f"fig_distance_sensitivity{suffix}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved fig_distance_sensitivity{suffix}.png")

    rank = sp.pivot(index="super_pop", columns="label", values="rank_by_mean")
    rank = rank.reindex(index=SUPERPOPS)
    fig, ax = plt.subplots(figsize=(9, 4))
    sns.heatmap(
        rank, cmap="viridis_r", annot=True, fmt=".0f", linewidths=0.5, ax=ax,
        cbar_kws={"label": "rank by mean score (1=highest)"},
    )
    ax.set_title(f"Population score-rank per PGS ({scope})", fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel("Super-population")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    fig.savefig(figs / f"fig_rank_instability{suffix}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved fig_rank_instability{suffix}.png")


if __name__ == "__main__":
    main()
