"""
Rebuild fine-map tiers with lead = max |z_meta| within LD block (not I2).

Avoids confounding fine_mapped with y_high_I2. Falls back to |beta_fixed|
then AF_max_diff if z_meta missing. SuSiE override still applied when present.
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
    p.add_argument("--susie-dir", default=str(ROOT / "data/labels/susie"))
    p.add_argument("--tag", default="genomewide_zlead")
    p.add_argument(
        "--out",
        default="",
        help="Default: data/labels/finemap_tiers_{tag}.parquet",
    )
    return p.parse_args()


def load_associated(master: Path) -> pd.DataFrame:
    import pyarrow.dataset as ds

    cols = [
        "variant_id",
        "trait",
        "I2",
        "y_high_I2",
        "associated",
        "sign_concordance",
        "z_meta",
        "beta_fixed",
        "AF_max_diff",
    ]
    dataset = ds.dataset(str(master), format="parquet")
    use = [c for c in cols if c in set(dataset.schema.names)]
    filt = ds.field("associated") == True  # noqa: E712
    # Write slim parquet once for reuse
    slim = ROOT / "data/labels/_tmp_associated_labels_z.parquet"
    if slim.exists() and slim.stat().st_size > 1_000_000:
        return pd.read_parquet(slim).drop_duplicates(["variant_id", "trait"])

    import pyarrow as pa
    import pyarrow.parquet as pq

    writer = None
    n = 0
    for batch in dataset.scanner(columns=use, filter=filt, batch_size=150_000).to_batches():
        t = pa.Table.from_batches([batch])
        if writer is None:
            writer = pq.ParquetWriter(str(slim), t.schema, compression="zstd")
        writer.write_table(t)
        n += t.num_rows
        if n % 1_000_000 < 150_000:
            print(f"  streamed {n:,}", flush=True)
    if writer:
        writer.close()
    print(f"  wrote {slim} ({n:,})", flush=True)
    return pd.read_parquet(slim).drop_duplicates(["variant_id", "trait"])


def lead_score(df: pd.DataFrame) -> pd.Series:
    if "z_meta" in df.columns:
        s = df["z_meta"].astype(float).abs()
        if s.notna().sum() > 1000:
            return s
    if "beta_fixed" in df.columns:
        s = df["beta_fixed"].astype(float).abs()
        if s.notna().sum() > 1000:
            return s
    return df.get("AF_max_diff", pd.Series(0, index=df.index)).astype(float).abs()


def heuristic_tiers(labels: pd.DataFrame, blocks: pd.DataFrame) -> pd.DataFrame:
    df = labels.merge(blocks[["variant_id", "ld_block"]], on="variant_id", how="left")
    df["ld_block"] = df["ld_block"].fillna("NA")
    df["lead_score"] = lead_score(df)
    df["rank_in_block"] = df.groupby(["trait", "ld_block"])["lead_score"].rank(
        ascending=False, method="first"
    )
    n_in_block = df.groupby(["trait", "ld_block"])["variant_id"].transform("count")
    df["finemap_tier"] = "tag_only"
    df.loc[(df["rank_in_block"] == 1) & (n_in_block >= 2), "finemap_tier"] = "fine_mapped"
    df.loc[n_in_block < 2, "finemap_tier"] = "ambiguous"
    df["tier_method"] = "ld_block_lead_abs_z_meta"
    return df[
        [
            "variant_id",
            "trait",
            "I2",
            "y_high_I2",
            "sign_concordance",
            "ld_block",
            "finemap_tier",
            "tier_method",
            "lead_score",
        ]
    ]


def apply_susie(tiers: pd.DataFrame, susie_dir: Path) -> pd.DataFrame:
    if not susie_dir.exists():
        return tiers
    files = list(susie_dir.glob("susie_*.parquet"))
    if not files:
        return tiers
    parts = []
    for f in files:
        name = f.stem.replace("susie_", "")
        trait = name.split("_")[0]
        s = pd.read_parquet(f)
        if "variant_id" not in s.columns:
            continue
        s["trait"] = trait
        s["in_cs"] = s.get("in_cs", s.get("pip", 0) > 0.5)
        if "fallback_mode" in s.columns:
            s["primary_cs"] = s["in_cs"] & (s["fallback_mode"].astype(str) == "signed_ld")
            s["pipeline_fallback"] = s["in_cs"] & (s["fallback_mode"].astype(str) != "signed_ld")
        else:
            # Legacy unsigned/identity outputs: never promote to primary fine_mapped
            s["primary_cs"] = False
            s["pipeline_fallback"] = s["in_cs"]
        parts.append(s[["variant_id", "trait", "primary_cs", "pipeline_fallback"]])
    if not parts:
        return tiers
    sus = pd.concat(parts, ignore_index=True)
    any_cs = sus.groupby(["variant_id", "trait"]).agg(
        primary_cs=("primary_cs", "any"),
        pipeline_fallback=("pipeline_fallback", "any"),
    ).reset_index()
    out = tiers.merge(any_cs, on=["variant_id", "trait"], how="left")
    out["primary_cs"] = out["primary_cs"].fillna(False)
    out["pipeline_fallback"] = out["pipeline_fallback"].fillna(False)
    out.loc[out["primary_cs"], "finemap_tier"] = "fine_mapped"
    out.loc[out["pipeline_fallback"] & ~out["primary_cs"], "finemap_tier"] = "pipeline_fallback"
    out.loc[
        ~out["primary_cs"] & ~out["pipeline_fallback"] & (out["finemap_tier"] != "ambiguous"),
        "finemap_tier",
    ] = "tag_only"
    out["tier_method"] = "susie_signed_ld_cs_only"
    return out.drop(columns=["primary_cs", "pipeline_fallback"], errors="ignore")


def main() -> None:
    args = parse_args()
    print("Loading associated labels (+ z_meta) ...", flush=True)
    labels = load_associated(Path(args.master))
    blocks = pd.read_parquet(args.ld_blocks, columns=["variant_id", "ld_block"])
    print(f"Labels {len(labels):,}  blocks {len(blocks):,}", flush=True)
    tiers = heuristic_tiers(labels, blocks)
    tiers = apply_susie(tiers, Path(args.susie_dir))
    out = Path(args.out) if args.out else ROOT / "data/labels" / f"finemap_tiers_{args.tag}.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    tiers.to_parquet(out, index=False)
    print(tiers.groupby(["trait", "finemap_tier"]).size().unstack(fill_value=0))
    # Sanity: fine_mapped should not be ~100% high-I2
    fm = tiers[tiers["finemap_tier"] == "fine_mapped"]
    if len(fm):
        print(
            f"fine_mapped pos_rate y_high_I2={fm['y_high_I2'].mean():.3f} "
            f"(n={len(fm):,}); tag_only={tiers.loc[tiers.finemap_tier=='tag_only','y_high_I2'].mean():.3f}"
        )
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
