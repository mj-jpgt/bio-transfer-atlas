"""
Phase B2: Non-additive / MOI proxy features from LD-block summaries.

Full i-LDSC can replace this when installed. Here we build a block-level
non-additivity proxy: residual variance of y_high_I2 after AF_LD_SEL-like
linear projection is unavailable offline, so we use LD entropy × AF_var
and within-block I2 dispersion as proxies, plus optional GenoBoost flags.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--master",
        default=str(ROOT / "data/modeling/master_variant_table_genomewide_genomewide.parquet"),
    )
    p.add_argument(
        "--ld-blocks",
        default=str(ROOT / "data/modeling/ld_block_assignments_genomewide.parquet"),
    )
    p.add_argument(
        "--genoboost",
        default=str(ROOT / "data/annotations/genoboost_moi.parquet"),
        help="Optional columns: variant_id, known_nonadditive (bool)",
    )
    p.add_argument(
        "--out",
        default=str(ROOT / "data/features/selection/nonadditive_proxy_features.parquet"),
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    import pyarrow.dataset as ds

    cols = [
        "variant_id",
        "trait",
        "I2",
        "associated",
        "AF_var",
        "LD_entropy",
        "FST_like",
        "LD_max_diff",
    ]
    # Prefer slim labels + join AF/LD from sample if full master too heavy
    slim = ROOT / "data/labels/_tmp_associated_labels.parquet"
    sample = ROOT / "data/modeling/_tmp_ldblock_associated_sample.parquet"
    if slim.exists() and sample.exists():
        import pyarrow.parquet as pq

        df = pd.read_parquet(slim, columns=["variant_id", "trait", "I2", "associated"])
        extra_cols = ["AF_var", "LD_entropy", "FST_like", "LD_max_diff"]
        names = set(pq.read_schema(sample).names)
        use_ex = ["variant_id"] + [c for c in extra_cols if c in names]
        extras = pd.read_parquet(sample, columns=use_ex).drop_duplicates("variant_id")
        df = df.merge(extras, on="variant_id", how="left")
    else:
        dataset = ds.dataset(str(args.master), format="parquet")
        use = [c for c in cols if c in set(dataset.schema.names)]
        filt = ds.field("associated") == True  # noqa: E712
        chunks = []
        for batch in dataset.scanner(columns=use, filter=filt, batch_size=200_000).to_batches():
            chunks.append(batch.to_pandas())
        df = pd.concat(chunks, ignore_index=True)
    blocks = pd.read_parquet(args.ld_blocks, columns=["variant_id", "ld_block"])
    df = df.merge(blocks, on="variant_id", how="left")

    # Per-variant proxies
    df["nonadd_proxy"] = (
        df.get("LD_entropy", 0).fillna(0).astype(float)
        * df.get("AF_var", 0).fillna(0).astype(float)
    )
    # Block-level I2 dispersion (heterogeneity of heterogeneity)
    g = df.groupby("ld_block")["I2"]
    df["i2_block_std"] = g.transform("std").fillna(0)
    df["i2_block_mean"] = g.transform("mean")
    df["nonadd_block_enrich"] = df["i2_block_std"] * df["nonadd_proxy"]

    out = (
        df.groupby("variant_id", as_index=False)
        .agg(
            nonadd_proxy=("nonadd_proxy", "mean"),
            i2_block_std=("i2_block_std", "mean"),
            nonadd_block_enrich=("nonadd_block_enrich", "mean"),
        )
    )
    gb = Path(args.genoboost)
    if gb.exists():
        gdf = pd.read_parquet(gb)
        out = out.merge(gdf, on="variant_id", how="left")
        if "known_nonadditive" in out.columns:
            out["known_nonadditive"] = out["known_nonadditive"].fillna(False).astype(bool)
    else:
        out["known_nonadditive"] = False

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out, index=False)
    print(f"Saved {args.out} ({len(out):,} variants)")
    print(out[["nonadd_proxy", "i2_block_std", "nonadd_block_enrich"]].describe())


if __name__ == "__main__":
    main()
