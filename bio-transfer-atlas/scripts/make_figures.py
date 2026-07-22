"""
Generate main atlas figures.

Usage:
    python scripts/make_figures.py \
      --metrics results/tables/population_shift_metrics.parquet \
      --out-dir results/figures
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import typer
from loguru import logger

app = typer.Typer(add_completion=False)


@app.command()
def main(
    metrics: Path = typer.Option(Path("results/tables/population_shift_metrics.parquet")),
    out_dir: Path = typer.Option(Path("results/figures")),
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(metrics)

    pivot = df.pivot_table(
        index="pgs_id", columns="population", values="standardized_shift", aggfunc="mean"
    )

    fig, ax = plt.subplots(figsize=(max(14, len(pivot.columns) * 0.4), max(6, len(pivot) * 0.5)))
    sns.heatmap(
        pivot,
        cmap="YlOrRd",
        linewidths=0.3,
        ax=ax,
        cbar_kws={"label": "Standardized shift"},
    )
    ax.set_title("Global PRS Score-Shift Atlas", fontsize=14, pad=12)
    ax.set_xlabel("Population")
    ax.set_ylabel("PGS ID")
    plt.tight_layout()

    out_path = out_dir / "atlas_heatmap.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.success(f"Atlas heatmap saved: {out_path}")


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="INFO")
    app()
