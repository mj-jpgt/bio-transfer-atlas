"""
Phase 18.3: Apply portability-risk interventions to PGS weights (genome-wide).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))

from intervention_common import (
    DUFFY_CHR,
    DUFFY_POS,
    DUFFY_WINDOW_BP,
    FILTER_FRACS,
    GW_INT_ROOT,
    GW_MASTER,
    GW_PREDS,
    GW_TAG,
    INTERVENTION_MODES,
    NEG_CONTROL_FRAC,
    PGS_IDS,
    PGS_TRAITS,
    REWEIGHT_EXP_K,
    ROOT,
    SEED,
    load_genomewide_weights,
    variant_chrom,
)

MDIR = ROOT / "data/modeling"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Apply intervention modes to genome-wide PGS weights.")
    p.add_argument("--predictions", default=str(GW_PREDS))
    p.add_argument("--master", default=str(GW_MASTER))
    p.add_argument("--pgs-root", default=str(ROOT / "data/processed/pgs_grch38"))
    p.add_argument("--out-root", default=str(GW_INT_ROOT))
    p.add_argument("--modes", default=",".join(INTERVENTION_MODES))
    p.add_argument("--seed", type=int, default=SEED)
    return p.parse_args()


def trait_risk(preds: pd.DataFrame, trait: str) -> pd.DataFrame:
    sub = preds[preds["trait"] == trait][["variant_id", "predicted_risk"]]
    return sub.drop_duplicates("variant_id", keep="first")


def join_weights_risk(
    weights: pd.DataFrame,
    preds: pd.DataFrame,
    trait: str,
    master_feats: pd.DataFrame | None = None,
) -> pd.DataFrame:
    risk = trait_risk(preds, trait)
    df = weights.merge(risk, on="variant_id", how="left")
    if master_feats is not None:
        feat_cols = [c for c in master_feats.columns if c not in df.columns]
        df = df.merge(master_feats[["variant_id"] + feat_cols], on="variant_id", how="left")
    df["predicted_risk"] = df["predicted_risk"].fillna(df["predicted_risk"].median())
    return df


def drop_top_frac(df: pd.DataFrame, col: str, frac: float, ascending: bool = False) -> pd.DataFrame:
    if len(df) == 0:
        return df
    n_drop = int(np.floor(len(df) * frac))
    if n_drop <= 0:
        return df.copy()
    order = df.sort_values(col, ascending=ascending)
    return order.iloc[n_drop:].copy()


def apply_mode(df: pd.DataFrame, mode: str, rng: np.random.Generator) -> pd.DataFrame:
    out = df.copy()
    if mode in FILTER_FRACS:
        out = drop_top_frac(out, "predicted_risk", FILTER_FRACS[mode], ascending=False)
    elif mode == "reweight_linear":
        out["effect_weight"] = out["effect_weight"] * (1.0 - out["predicted_risk"])
    elif mode == "reweight_exp":
        out["effect_weight"] = out["effect_weight"] * np.exp(-REWEIGHT_EXP_K * out["predicted_risk"])
    elif mode == "flag":
        pass
    elif mode == "random":
        out = out.sample(frac=1.0 - NEG_CONTROL_FRAC, random_state=int(rng.integers(0, 2**31 - 1)))
    elif mode == "fst":
        out = drop_top_frac(out, "FST_like", NEG_CONTROL_FRAC, ascending=False)
    elif mode == "ld":
        out = drop_top_frac(out, "LD_EUR", NEG_CONTROL_FRAC, ascending=False)
    elif mode == "maf":
        maf_col = "AF_EUR" if "AF_EUR" in out.columns else "AF_max_diff"
        out = drop_top_frac(out, maf_col, NEG_CONTROL_FRAC, ascending=True)
    elif mode == "duffy_gate":
        # Drop variants in ACKR1/Duffy-null window (chr1 ±200kb around rs2814778)
        def in_duffy(vid: str) -> bool:
            parts = str(vid).split(":")
            if len(parts) < 2:
                return False
            try:
                chrom = variant_chrom(vid)
                pos = int(float(parts[1]))
            except Exception:
                return False
            return chrom == DUFFY_CHR and abs(pos - DUFFY_POS) <= DUFFY_WINDOW_BP

        out = out[~out["variant_id"].map(in_duffy)].copy()
    else:
        raise ValueError(f"Unknown mode: {mode}")
    return out


def weight_retained_frac(before: pd.DataFrame, after: pd.DataFrame) -> float:
    denom = before["effect_weight"].abs().sum()
    if denom == 0:
        return float("nan")
    return float(after["effect_weight"].abs().sum() / denom)


def write_plink_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    out = df[["variant_id", "effect_allele", "effect_weight"]].copy()
    out.to_csv(path, sep="\t", index=False)


def load_master_feats(master_path: Path) -> pd.DataFrame:
    import pyarrow.dataset as ds

    cols = ["variant_id", "FST_like", "LD_EUR", "AF_EUR"]
    dataset = ds.dataset(str(master_path), format="parquet")
    scanner = dataset.scanner(columns=cols, batch_size=1_000_000)
    chunks = []
    for batch in scanner.to_batches():
        chunks.append(batch.to_pandas())
    if not chunks:
        return pd.DataFrame(columns=cols)
    # Master is variant×trait; keep one feature row per variant
    return pd.concat(chunks, ignore_index=True).drop_duplicates("variant_id", keep="first")


def main() -> None:
    args = parse_args()
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    out_root = Path(args.out_root)
    pgs_root = Path(args.pgs_root)

    preds = pd.read_parquet(args.predictions)
    print("Loading master features (deduped by variant_id) ...", flush=True)
    master_feats = load_master_feats(Path(args.master))
    print(f"  {len(master_feats):,} unique variants with AF/LD/FST", flush=True)

    rng = np.random.default_rng(args.seed)
    manifest_rows = []

    for pgs_id in PGS_IDS:
        trait = PGS_TRAITS[pgs_id]
        weights = load_genomewide_weights(pgs_id, pgs_root)
        weights = weights.drop_duplicates("variant_id", keep="first")
        n_in = len(weights)
        base = join_weights_risk(weights, preds, trait, master_feats)
        base = base.drop_duplicates("variant_id", keep="first")
        print(f"{pgs_id}: {n_in:,} weights after dedup", flush=True)

        for mode in modes:
            modified = apply_mode(base, mode, rng)
            n_out = len(modified)
            out_path = out_root / pgs_id / f"{mode}.tsv"
            write_plink_tsv(modified, out_path)
            manifest_rows.append({
                "pgs_id": pgs_id,
                "trait": trait,
                "mode": mode,
                "n_variants_in": n_in,
                "n_variants_out": n_out,
                "frac_retained": n_out / n_in if n_in else float("nan"),
                "weight_retained_frac": weight_retained_frac(weights, modified),
                "out_path": str(out_path),
            })
            print(f"  {pgs_id} {mode}: {n_in} -> {n_out} variants")

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = out_root / f"intervention_manifest.{GW_TAG}.parquet"
    manifest.to_parquet(manifest_path, index=False)
    meta = {"modes": modes, "seed": args.seed, "predictions": args.predictions, "scope": GW_TAG}
    (out_root / f"intervention_manifest.{GW_TAG}.json").write_text(json.dumps(meta, indent=2))
    print(f"\nManifest -> {manifest_path}")


if __name__ == "__main__":
    main()
