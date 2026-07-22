"""
Compute population-level instability metrics from score matrix + sample metadata.

Usage:
    python scripts/compute_metrics.py \
      --score-matrix data/processed/scores/PGS000XXX/all_chrom.sscore \
      --sample-meta  data/processed/ancestry/sample_metadata.parquet \
      --pgs-id       PGS000XXX \
      --out          data/processed/metrics/PGS000XXX/population_shift.parquet \
      --n-boot       1000 \
      --seed         719
"""

import sys
from pathlib import Path

import pandas as pd
import typer
from loguru import logger

from bta.metrics.shifts import standardized_shift, bootstrap_shifts

app = typer.Typer(add_completion=False)


@app.command()
def main(
    score_matrix: Path = typer.Option(...),
    sample_meta: Path = typer.Option(...),
    pgs_id: str = typer.Option(...),
    out: Path = typer.Option(...),
    n_boot: int = typer.Option(1000),
    seed: int = typer.Option(719),
    reference: str = typer.Option("EUR", help="Reference super-population"),
) -> None:
    scores_df = pd.read_csv(score_matrix, sep="\t")
    meta = pd.read_parquet(sample_meta)

    score_col = next(
        (c for c in scores_df.columns if "score" in c.lower() or "sum" in c.lower()),
        scores_df.columns[-1],
    )
    merged = scores_df.merge(meta[["sample", "population", "super_population"]], on="sample", how="inner")

    scores = merged[score_col]
    labels = merged["population"]
    super_labels = merged["super_population"]

    pop_shifts = standardized_shift(scores, labels, reference=None)
    super_shifts = standardized_shift(scores, super_labels, reference=reference)

    boot_ci = bootstrap_shifts(scores, labels, reference=None, n_boot=n_boot, seed=seed)

    combined = pop_shifts.merge(boot_ci, on="population", how="left")
    combined["pgs_id"] = pgs_id

    out.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(out, index=False)
    logger.success(f"Metrics saved: {out}  ({len(combined)} populations)")


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="INFO")
    app()
