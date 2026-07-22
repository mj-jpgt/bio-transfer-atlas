"""
Phase B1: Fine-mapping tier labels (SuSiE output if present; else LD-block lead heuristic).

Tier:
  fine_mapped — lead variant in LD block (or PIP>0.5 / in CS from SuSiE)
  tag_only — associated but never lead/CS
  ambiguous — not associated or insufficient block mates
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SEED = 719


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
        "--susie-dir",
        default=str(ROOT / "data/labels/susie"),
        help="Optional dir with susie_{trait}_{anc}.parquet columns: variant_id,pip,in_cs",
    )
    p.add_argument("--tag", default="genomewide")
    return p.parse_args()


def load_labels(master: Path) -> pd.DataFrame:
    # Prefer slim on-disk associated labels if present (low RAM)
    slim = ROOT / "data/labels/_tmp_associated_labels.parquet"
    if slim.exists():
        df = pd.read_parquet(slim)
        return df.drop_duplicates(["variant_id", "trait"])
    import pyarrow.dataset as ds

    cols = ["variant_id", "trait", "I2", "y_high_I2", "associated", "sign_concordance"]
    dataset = ds.dataset(str(master), format="parquet")
    use = [c for c in cols if c in set(dataset.schema.names)]
    filt = ds.field("associated") == True  # noqa: E712
    chunks = []
    for batch in dataset.scanner(columns=use, filter=filt, batch_size=200_000).to_batches():
        chunks.append(batch.to_pandas())
    return pd.concat(chunks, ignore_index=True).drop_duplicates(["variant_id", "trait"])


def heuristic_tiers(labels: pd.DataFrame, blocks: pd.DataFrame) -> pd.DataFrame:
    df = labels.merge(blocks[["variant_id", "ld_block"]], on="variant_id", how="left")
    df["ld_block"] = df["ld_block"].fillna("NA")
    # Lead = highest I2 within block×trait among associated (proxy for causal heterogeneity locus)
    df["rank_in_block"] = df.groupby(["trait", "ld_block"])["I2"].rank(
        ascending=False, method="first"
    )
    n_in_block = df.groupby(["trait", "ld_block"])["variant_id"].transform("count")
    df["finemap_tier"] = "tag_only"
    df.loc[(df["rank_in_block"] == 1) & (n_in_block >= 2), "finemap_tier"] = "fine_mapped"
    df.loc[n_in_block < 2, "finemap_tier"] = "ambiguous"
    df["tier_method"] = "ld_block_lead_I2"
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
        # Primary fine_mapped requires signed LD SuSiE only
        if "fallback_mode" in s.columns:
            primary = s.query("fallback_mode == 'signed_ld'")
            s["primary_cs"] = False
            s.loc[primary.index, "primary_cs"] = primary["in_cs"].astype(bool)
            s["pipeline_fallback"] = s["in_cs"] & (s["fallback_mode"].astype(str) != "signed_ld")
        else:
            # Legacy files without fallback_mode: treat as non-primary
            s["primary_cs"] = False
            s["pipeline_fallback"] = s["in_cs"]
        parts.append(s[["variant_id", "trait", "primary_cs", "pipeline_fallback"]])
    if not parts:
        return tiers
    sus = pd.concat(parts, ignore_index=True)
    any_cs = sus.groupby(["variant_id", "trait"], as_index=False).agg(
        primary_cs=("primary_cs", "any"),
        pipeline_fallback=("pipeline_fallback", "any"),
    )
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
    print("Loading labels ...", flush=True)
    labels = load_labels(Path(args.master))
    blocks = pd.read_parquet(args.ld_blocks, columns=["variant_id", "ld_block"])
    print(f"Labels {len(labels):,}  blocks {len(blocks):,}", flush=True)
    tiers = heuristic_tiers(labels, blocks)
    tiers = apply_susie(tiers, Path(args.susie_dir))
    out = ROOT / "data/labels" / f"finemap_tiers_{args.tag}.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    tiers.to_parquet(out, index=False)
    print(tiers.groupby(["trait", "finemap_tier"]).size().unstack(fill_value=0))
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
