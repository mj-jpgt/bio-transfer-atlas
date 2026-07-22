"""
FAIRGEN-Open Stage 10: Master variant table + modeling splits (chr22)
======================================================================
Joins cross-ancestry concordance labels with AF / LD / selection features,
defines feature groups, and assigns leakage-safe splits.

Splits:
  split_trait   : biological-domain generalization
                    train = CAD + T2D, val = BMI, test = LDL
  split_variant : rigorous group split by variant_id (no variant in >1 fold)
                    70/15/15 train/val/test

Outputs:
  data/modeling/master_variant_table.parquet
  data/modeling/feature_groups.json
"""
from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FEAT = ROOT / "data/features"
OUT_DIR = ROOT / "data/modeling"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SEED = 719


def normalize_chr(value: str) -> str:
    s = str(value).strip()
    return s[3:] if s.lower().startswith("chr") else s


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22", help="Chromosome number/name (e.g. 22 or chr22)")
    return p.parse_args()


def maybe_write_legacy(chrom: str, src: Path, legacy_name: str) -> None:
    if chrom == "22":
        pd.read_parquet(src).to_parquet(OUT_DIR / legacy_name, index=False)


def resolve_label_path(chrom: str) -> Path:
    cands = [
        ROOT / "data/labels" / f"gwas_concordance_labels_multisource.chr{chrom}.parquet",
        ROOT / "data/labels" / f"gwas_concordance_labels.chr{chrom}.parquet",
    ]
    if chrom == "22":
        cands.extend(
            [
                ROOT / "data/labels/gwas_concordance_labels_multisource.parquet",
                ROOT / "data/labels/gwas_concordance_labels.parquet",
            ]
        )
    for p in cands:
        if p.exists():
            return p
    raise FileNotFoundError(f"No label parquet found for chr{chrom}")


def resolve_feature_path(kind: str, chrom: str, legacy: str) -> Path:
    p = FEAT / kind / f"{legacy}.chr{chrom}.parquet"
    if p.exists():
        return p
    if chrom == "22":
        legacy_path = FEAT / kind / f"{legacy}.parquet"
        if legacy_path.exists():
            return legacy_path
    raise FileNotFoundError(f"Missing {kind} features for chr{chrom}: {p}")


def resolve_optional_selection_feature(stem: str, chrom: str) -> Path | None:
    p = FEAT / "selection" / f"{stem}.chr{chrom}.parquet"
    if p.exists():
        return p
    if chrom == "22":
        legacy = FEAT / "selection" / f"{stem}.parquet"
        if legacy.exists():
            return legacy
    return None


def main() -> None:
    args = parse_args()
    chrom = normalize_chr(args.chrom)
    rng = np.random.default_rng(SEED)

    labels = resolve_label_path(chrom)
    af_path = resolve_feature_path("af", chrom, "af_divergence_features")
    ld_path = resolve_feature_path("ld", chrom, "ld_divergence_features")
    sel_path = resolve_feature_path("selection", chrom, "selection_turnover_features")

    lab = pd.read_parquet(labels)
    af = pd.read_parquet(af_path)
    ld = pd.read_parquet(ld_path)
    sel = pd.read_parquet(sel_path)
    constraint_path = resolve_optional_selection_feature("constraint_features", chrom)
    conservation_path = resolve_optional_selection_feature("conservation_features", chrom)
    constraint = pd.read_parquet(constraint_path) if constraint_path is not None else None
    conservation = pd.read_parquet(conservation_path) if conservation_path is not None else None
    print(
        f"chr{chrom} labels={len(lab):,} af={len(af):,} ld={len(ld):,} sel={len(sel):,} "
        f"constraint={'yes' if constraint is not None else 'no'} "
        f"conservation={'yes' if conservation is not None else 'no'}"
    )

    af_drop = [c for c in ["REF", "ALT", "AF_global"] if c in af.columns]
    master = (
        lab.merge(af.drop(columns=af_drop), on="variant_id", how="inner", validate="many_to_one")
        .merge(ld, on="variant_id", how="inner", validate="many_to_one")
        .merge(sel, on="variant_id", how="inner", validate="many_to_one")
    )
    if constraint is not None:
        master = master.merge(constraint, on="variant_id", how="left", validate="many_to_one")
    if conservation is not None:
        master = master.merge(conservation, on="variant_id", how="left", validate="many_to_one")
    print(f"master rows after merge: {len(master):,}")

    af_feats = [c for c in af.columns if c.startswith(("AF_", "FST_like", "AIS")) and c != "AF_global"]
    ld_feats = [c for c in ld.columns if c.startswith("LD_")]
    sel_feats = [c for c in sel.columns if c.startswith(("FST_", "PBS_", "DAF_")) and c != "PBS_argmax_pop"]
    has_constraint = constraint is not None
    has_conservation = conservation is not None
    del lab, af, ld, sel, constraint, conservation
    gc.collect()
    if has_constraint:
        sel_feats += [c for c in ["LOEUF", "pLI", "mis_z"] if c in master.columns]
    if has_conservation:
        sel_feats += [c for c in ["phyloP", "phyloP_win"] if c in master.columns]
    sel_feats = sorted(set(sel_feats))
    fst_only = ["FST_like"] if "FST_like" in master.columns else []

    all_feats = sorted(set(af_feats + ld_feats + sel_feats))
    for c in all_feats:
        master[c] = pd.to_numeric(master[c], errors="coerce").replace([np.inf, -np.inf], np.nan)

    feature_groups = {
        "FST": fst_only,
        "AF": af_feats,
        "LD": ld_feats,
        "SEL": sel_feats,
        "AF_LD": af_feats + ld_feats,
        "AF_SEL": af_feats + sel_feats,
        "LD_SEL": ld_feats + sel_feats,
        "AF_LD_SEL": af_feats + ld_feats + sel_feats,
    }
    fg_path = OUT_DIR / f"feature_groups.chr{chrom}.json"
    fg_path.write_text(json.dumps(feature_groups, indent=2), encoding="utf-8")
    if chrom == "22":
        (OUT_DIR / "feature_groups.json").write_text(json.dumps(feature_groups, indent=2), encoding="utf-8")
    print("Feature group sizes:", {k: len(v) for k, v in feature_groups.items()})

    master["y_sign_discordant"] = (master["sign_concordance"] < 1.0).astype(int)
    master["y_high_I2"] = (master["I2"] > 0.25).astype(int)
    master["y_risk_class"] = master["risk_class"].astype(int)

    trait_map = {"CAD": "train", "T2D": "train", "BMI": "val", "LDL": "test"}
    master["split_trait"] = master["trait"].map(trait_map)

    uniq = master["variant_id"].drop_duplicates().to_numpy()
    perm = rng.permutation(len(uniq))
    n_tr = int(0.70 * len(uniq))
    n_va = int(0.15 * len(uniq))
    folds = np.full(len(uniq), "test", dtype=object)
    folds[perm[:n_tr]] = "train"
    folds[perm[n_tr : n_tr + n_va]] = "val"
    split_df = pd.DataFrame({"variant_id": uniq, "split_variant": folds})
    split_df["split_variant"] = split_df["split_variant"].astype("category")
    master = master.merge(split_df, on="variant_id", how="left")
    del split_df, uniq, perm, folds
    gc.collect()

    # Downcast dtypes and write via ParquetWriter in chunks to avoid ArrowMemoryError
    import pyarrow as pa
    import pyarrow.parquet as pq
    for col in master.select_dtypes(include="float64").columns:
        master[col] = master[col].astype("float32")
    for col in ("trait", "source_group", "risk_class", "split_trait", "split_variant"):
        if col in master.columns:
            master[col] = master[col].astype("category")
    out = OUT_DIR / f"master_variant_table.chr{chrom}.parquet"
    chunk_size = 500_000
    writer = None
    for start in range(0, len(master), chunk_size):
        chunk = master.iloc[start : start + chunk_size]
        table = pa.Table.from_pandas(chunk, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(str(out), table.schema)
        writer.write_table(table)
    if writer:
        writer.close()
    maybe_write_legacy(chrom, out, "master_variant_table.parquet")
    print(f"\nSaved {out} ({len(master):,} rows × {master.shape[1]} cols)")

    print("\nTarget prevalence (all rows):")
    print(f"  y_sign_discordant: {master['y_sign_discordant'].mean():.3f}")
    print(f"  y_high_I2        : {master['y_high_I2'].mean():.3f}")
    print("  risk_class:", master["y_risk_class"].value_counts(normalize=True).round(3).to_dict())

    assoc = master[master["associated"]]
    print(f"\nAssociated subset: {len(assoc):,} rows")
    print(f"  y_sign_discordant: {assoc['y_sign_discordant'].mean():.3f}")
    print(f"  y_high_I2        : {assoc['y_high_I2'].mean():.3f}")
    print("  risk_class:", assoc["y_risk_class"].value_counts(normalize=True).round(3).to_dict())
    print("\nsplit_trait counts:", master["split_trait"].value_counts().to_dict())


if __name__ == "__main__":
    main()
