"""
Merge per-chromosome PLINK2 .sscore files into one all-chromosome score matrix.

Usage:
    python scripts/calculate_scores.py \
      --pgs-id PGS000XXX \
      --chroms 22 \
      --score-dir data/processed/scores \
      --out data/processed/scores/PGS000XXX/all_chrom.sscore
"""

import sys
from pathlib import Path

import pandas as pd
import typer
from loguru import logger

app = typer.Typer(add_completion=False)


@app.command()
def main(
    pgs_id: str = typer.Option(...),
    chroms: str = typer.Option("22", help="Comma-separated chromosome list"),
    score_dir: Path = typer.Option(Path("data/processed/scores")),
    out: Path = typer.Option(...),
) -> None:
    chrom_list = [c.strip() for c in chroms.split(",")]
    frames = []
    for chrom in chrom_list:
        p = score_dir / pgs_id / f"chr{chrom}.sscore"
        if not p.exists():
            logger.warning(f"Missing: {p}")
            continue
        df = pd.read_csv(p, sep="\t")
        df.columns = [c.lower() for c in df.columns]
        df["chrom"] = chrom
        frames.append(df)

    if not frames:
        logger.error("No sscore files found.")
        raise SystemExit(1)

    merged = pd.concat(frames, ignore_index=True)

    id_col = next((c for c in merged.columns if "iid" in c or "fid" in c), None)
    if id_col:
        merged = merged.rename(columns={id_col: "sample"})

    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, sep="\t", index=False)
    logger.success(f"Score matrix: {out}  ({len(merged)} rows)")


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="INFO")
    app()
