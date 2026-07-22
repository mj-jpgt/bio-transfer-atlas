"""
Merge 1000G sample metadata with PLINK2 PCA eigenvec output.
Produces:
  data/processed/ancestry/sample_metadata.parquet
  data/processed/ancestry/ancestry_pcs.parquet

Usage (called by Snakemake or directly):
    python scripts/compute_pcs.py \
      --panel  data/raw/1000g/metadata/integrated_call_samples_v3.20130502.ALL.panel \
      --eigenvec data/processed/1000g/chr22_pca.eigenvec \
      --out-meta data/processed/ancestry/sample_metadata.parquet \
      --out-pcs  data/processed/ancestry/ancestry_pcs.parquet
"""

import sys
from pathlib import Path

import pandas as pd
import typer
from loguru import logger

app = typer.Typer(add_completion=False)


@app.command()
def main(
    panel: Path = typer.Option(...),
    eigenvec: Path = typer.Option(...),
    out_meta: Path = typer.Option(...),
    out_pcs: Path = typer.Option(...),
) -> None:
    meta = pd.read_csv(panel, sep="\t")
    meta.columns = [c.lower() for c in meta.columns]
    meta = meta.rename(columns={"sample": "sample", "pop": "population", "super_pop": "super_population", "gender": "gender"})

    pcs_raw = pd.read_csv(eigenvec, sep="\t")
    pcs_raw.columns = [c.lower() for c in pcs_raw.columns]

    id_col = "#iid" if "#iid" in pcs_raw.columns else "iid"
    pcs_raw = pcs_raw.rename(columns={id_col: "sample"})
    if "fid" in pcs_raw.columns:
        pcs_raw = pcs_raw.drop(columns=["fid"])

    pc_cols = [c for c in pcs_raw.columns if c.startswith("pc")]
    pcs = pcs_raw[["sample"] + pc_cols]

    merged_meta = meta.merge(pcs[["sample"]], on="sample", how="inner")
    pcs_merged = pcs.merge(meta[["sample", "population", "super_population"]], on="sample", how="left")

    out_meta.parent.mkdir(parents=True, exist_ok=True)
    out_pcs.parent.mkdir(parents=True, exist_ok=True)

    merged_meta.to_parquet(out_meta, index=False)
    pcs_merged.to_parquet(out_pcs, index=False)

    logger.success(f"sample_metadata: {out_meta}  ({len(merged_meta)} samples)")
    logger.success(f"ancestry_pcs:    {out_pcs}  ({len(pcs_merged)} samples, {len(pc_cols)} PCs)")


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="INFO")
    app()
