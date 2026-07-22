"""
Phase A2: Assign LD-block holdout splits (1 Mb genomic bins; LDetect if present).

Writes:
  data/modeling/ld_block_assignments_genomewide.parquet
  and updates a sidecar split table for evaluation.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SEED = 719
RNG = np.random.default_rng(SEED)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Assign LD-block train/test splits.")
    p.add_argument(
        "--master",
        default=str(ROOT / "data/modeling/master_variant_table_genomewide_genomewide.parquet"),
    )
    p.add_argument("--block-mb", type=float, default=1.0)
    p.add_argument("--test-frac", type=float, default=0.2)
    p.add_argument(
        "--ldetect",
        default=str(ROOT / "data/annotations/ldetect_grch38.bed"),
        help="Optional BED: chrom start end (0-based)",
    )
    p.add_argument(
        "--out",
        default=str(ROOT / "data/modeling/ld_block_assignments_genomewide.parquet"),
    )
    return p.parse_args()


def parse_vid(vid: str) -> tuple[str, int] | None:
    parts = str(vid).split(":")
    if len(parts) < 2:
        return None
    try:
        return parts[0], int(parts[1])
    except ValueError:
        return None


def load_unique_variants(master: Path) -> pd.DataFrame:
    import pyarrow.dataset as ds

    dataset = ds.dataset(str(master), format="parquet")
    filt = ds.field("associated") == True  # noqa: E712
    chunks = []
    for batch in dataset.scanner(columns=["variant_id"], filter=filt, batch_size=1_000_000).to_batches():
        chunks.append(batch.to_pandas())
    vids = pd.concat(chunks, ignore_index=True)["variant_id"].drop_duplicates()
    rows = []
    for vid in vids:
        parsed = parse_vid(vid)
        if parsed is None:
            continue
        chrom, pos = parsed
        rows.append({"variant_id": vid, "chrom": chrom, "pos": pos})
    return pd.DataFrame(rows)


def assign_mb_blocks(df: pd.DataFrame, block_mb: float) -> pd.Series:
    bp = int(block_mb * 1_000_000)
    return df["chrom"].astype(str) + ":" + (df["pos"] // bp).astype(str)


def assign_ldetect(df: pd.DataFrame, bed: Path) -> pd.Series | None:
    if not bed.exists():
        return None
    blocks = pd.read_csv(
        bed, sep="\t", header=None, names=["chrom", "start", "end"], comment="#"
    )
    blocks["chrom"] = blocks["chrom"].astype(str).str.replace("^chr", "", regex=True)
    # interval join per chrom
    out = pd.Series(index=df.index, dtype=object)
    for chrom, sub in df.groupby("chrom"):
        b = blocks[blocks["chrom"] == str(chrom)].sort_values("start")
        if b.empty:
            out.loc[sub.index] = str(chrom) + ":NA"
            continue
        starts = b["start"].to_numpy()
        ends = b["end"].to_numpy()
        labels = (b["chrom"].astype(str) + ":" + b["start"].astype(str) + "-" + b["end"].astype(str)).to_numpy()
        pos = sub["pos"].to_numpy()
        idx = np.searchsorted(starts, pos, side="right") - 1
        assigned = []
        for i, p in enumerate(pos):
            j = idx[i]
            if j >= 0 and j < len(ends) and starts[j] <= p < ends[j]:
                assigned.append(labels[j])
            else:
                assigned.append(f"{chrom}:unassigned:{p // 1_000_000}")
        out.loc[sub.index] = assigned
    return out


def main() -> None:
    args = parse_args()
    print("Loading unique associated variants ...", flush=True)
    df = load_unique_variants(Path(args.master))
    print(f"  {len(df):,} unique variants", flush=True)

    ldetect = Path(args.ldetect)
    block_ids = assign_ldetect(df, ldetect)
    method = "ldetect"
    if block_ids is None:
        print(f"LDetect BED not found at {ldetect}; using {args.block_mb} Mb bins", flush=True)
        block_ids = assign_mb_blocks(df, args.block_mb)
        method = f"mb_bin_{args.block_mb}"
    df["ld_block"] = block_ids

    blocks = df["ld_block"].drop_duplicates().to_numpy()
    RNG.shuffle(blocks)
    n_test = max(1, int(len(blocks) * args.test_frac))
    test_blocks = set(blocks[:n_test])
    df["split_ld_block"] = df["ld_block"].map(lambda b: "test" if b in test_blocks else "train")

    # Sanity: no block in both
    tr_b = set(df.loc[df["split_ld_block"] == "train", "ld_block"])
    te_b = set(df.loc[df["split_ld_block"] == "test", "ld_block"])
    assert tr_b.isdisjoint(te_b), "LD-block split leakage"

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df[["variant_id", "chrom", "pos", "ld_block", "split_ld_block"]].to_parquet(out, index=False)
    meta = {
        "method": method,
        "n_variants": len(df),
        "n_blocks": int(df["ld_block"].nunique()),
        "n_train_blocks": len(tr_b),
        "n_test_blocks": len(te_b),
        "test_frac_variants": float((df["split_ld_block"] == "test").mean()),
    }
    meta_path = out.with_suffix(".json")
    import json

    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
