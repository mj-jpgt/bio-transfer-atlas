"""
Phase 18.6: Per-individual PGS reliability flag (genome-wide, multi-chromosome).

Output is a descriptive reliability flag only — not a clinically calibrated score.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from intervention_common import (
    GW_PREDS,
    GW_TAG,
    PGS_IDS,
    PGS_TRAITS,
    PLINK2,
    ROOT,
    available_score_chroms,
    load_genomewide_weights,
    pfile_for_chrom,
    variant_chrom,
)

OUT_DEFAULT = ROOT / f"results/tables/per_patient_confidence.{GW_TAG}.parquet"
BATCH_SIZE = 500
PLINK_MEMORY_MB = 640


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Per-patient intervention confidence flags (GW).")
    p.add_argument("--predictions", default=str(GW_PREDS))
    p.add_argument("--chroms", default="", help="Comma list; default available score pgens")
    p.add_argument("--out", default=str(OUT_DEFAULT))
    p.add_argument("--decile", type=float, default=0.10)
    p.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    p.add_argument(
        "--max-variants",
        type=int,
        default=20_000,
        help="Cap variants by |effect_weight| for runtime (descriptive flag only)",
    )
    return p.parse_args()


def parse_chrom_list(spec: str) -> list[str]:
    if not spec.strip():
        return available_score_chroms()
    return [str(int(x.strip())) for x in spec.split(",") if x.strip()]


def trait_risk_map(preds: pd.DataFrame, trait: str) -> dict[str, float]:
    sub = preds[preds["trait"] == trait][["variant_id", "predicted_risk"]]
    sub = sub.drop_duplicates("variant_id", keep="first")
    return dict(zip(sub["variant_id"], 1.0 - sub["predicted_risk"]))


def plink_raw_variant_id(col: str) -> str:
    return col.rsplit("_", 1)[0]


def export_dosage_batch(
    pfile: Path,
    variant_ids: list[str],
    work_dir: Path,
    tag: str,
) -> tuple[np.ndarray, list[str], np.ndarray, list[str]]:
    extract = work_dir / f"extract_{tag}.txt"
    extract.write_text("\n".join(variant_ids) + "\n")
    out_prefix = work_dir / f"dosages_{tag}"
    cmd = [
        PLINK2,
        "--pfile", str(pfile),
        "--memory", str(PLINK_MEMORY_MB),
        "--extract", str(extract),
        "--export", "A", "ref-first",
        "--out", str(out_prefix),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"plink2 export failed ({tag}):\n{r.stderr[-500:]}")
    raw = pd.read_csv(Path(str(out_prefix) + ".raw"), sep="\t")
    sample_ids = raw["IID"].astype(str).to_numpy()
    meta_cols = {"FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"}
    raw_cols = [c for c in raw.columns if c not in meta_cols]
    mapped_ids = [plink_raw_variant_id(c) for c in raw_cols]
    dose = raw[raw_cols].to_numpy(dtype=np.float32)
    return dose, mapped_ids, sample_ids, raw_cols


def update_top_decile(
    top_contrib: np.ndarray | None,
    top_conf: np.ndarray | None,
    batch_contrib: np.ndarray,
    batch_conf: np.ndarray,
    n_top: int,
) -> tuple[np.ndarray, np.ndarray]:
    if top_contrib is None:
        top_contrib = batch_contrib
        top_conf = batch_conf
    else:
        top_contrib = np.concatenate([top_contrib, batch_contrib], axis=1)
        top_conf = np.concatenate([top_conf, batch_conf], axis=1)

    if top_contrib.shape[1] > n_top:
        idx = np.argpartition(-top_contrib, n_top - 1, axis=1)[:, :n_top]
        top_contrib = np.take_along_axis(top_contrib, idx, axis=1)
        top_conf = np.take_along_axis(top_conf, idx, axis=1)
    return top_contrib, top_conf


def filter_weights_for_chroms(weights: pd.DataFrame, chroms: list[str], max_variants: int) -> pd.DataFrame:
    chrom_set = set(chroms)
    w = weights.dropna(subset=["effect_weight"]).copy()
    w["chrom"] = w["variant_id"].map(variant_chrom)
    w = w[w["chrom"].isin(chrom_set)]
    if max_variants > 0 and len(w) > max_variants:
        w = w.assign(_abs=w["effect_weight"].abs()).nlargest(max_variants, "_abs").drop(columns=["_abs"])
    return w.drop(columns=["chrom"], errors="ignore")


def per_patient_for_pgs(
    chroms: list[str],
    weights: pd.DataFrame,
    conf_map: dict[str, float],
    decile: float,
    batch_size: int,
    work_dir: Path,
    pgs_id: str,
) -> pd.DataFrame:
    w = weights.dropna(subset=["effect_weight"]).copy()
    variants = w["variant_id"].tolist()
    if not variants:
        return pd.DataFrame(columns=["sample_id", "score_confidence"])

    wt_series = w.set_index("variant_id")["effect_weight"]
    sample_ids: np.ndarray | None = None
    n_top = max(1, int(np.ceil(len(variants) * decile)))
    top_contrib: np.ndarray | None = None
    top_conf: np.ndarray | None = None
    batch_idx = 0

    for chrom in chroms:
        pfile = pfile_for_chrom(chrom)
        chrom_vars = [v for v in variants if variant_chrom(v) == chrom]
        if not chrom_vars:
            continue
        for i in range(0, len(chrom_vars), batch_size):
            batch_vars = chrom_vars[i : i + batch_size]
            dose, mapped_ids, sids, _ = export_dosage_batch(
                pfile, batch_vars, work_dir, f"{pgs_id}_{chrom}_{batch_idx}"
            )
            batch_idx += 1
            if sample_ids is None:
                sample_ids = sids
            elif not np.array_equal(sample_ids, sids):
                raise RuntimeError(f"Sample order mismatch in {pgs_id} chr{chrom}")

            wt = wt_series.reindex(mapped_ids).to_numpy(dtype=np.float32)
            conf = np.array([conf_map.get(v, np.nan) for v in mapped_ids], dtype=np.float32)
            batch_contrib = np.abs(dose * wt[np.newaxis, :])
            batch_conf = np.broadcast_to(conf[np.newaxis, :], batch_contrib.shape)
            top_contrib, top_conf = update_top_decile(
                top_contrib, top_conf, batch_contrib, batch_conf, n_top
            )

    if sample_ids is None or top_contrib is None:
        return pd.DataFrame(columns=["sample_id", "score_confidence"])
    score_conf = np.nanmean(top_conf, axis=1)
    return pd.DataFrame({"sample_id": sample_ids, "score_confidence": score_conf})


def main() -> None:
    args = parse_args()
    chroms = parse_chrom_list(args.chroms)
    preds = pd.read_parquet(args.predictions)
    panel = pd.read_parquet(ROOT / "data/processed/sample_metadata_grch38.parquet")

    rows = []
    with tempfile.TemporaryDirectory(prefix="ppc_gw_") as tmp:
        work = Path(tmp)
        for pgs_id in PGS_IDS:
            trait = PGS_TRAITS[pgs_id]
            print(f"{pgs_id} ({trait}) ...")
            weights = load_genomewide_weights(pgs_id)
            weights = weights.drop_duplicates("variant_id", keep="first")
            weights = filter_weights_for_chroms(weights, chroms, args.max_variants)
            print(f"  {len(weights):,} variants on chroms {chroms} (capped)")
            conf_map = trait_risk_map(preds, trait)
            ppc = per_patient_for_pgs(
                chroms, weights, conf_map, args.decile, args.batch_size, work, pgs_id
            )
            ppc["pgs_id"] = pgs_id
            rows.append(ppc)
            if len(ppc):
                print(f"  mean confidence = {ppc['score_confidence'].mean():.4f}")

    out = pd.concat(rows, ignore_index=True)
    out = out.merge(panel[["sample_id", "super_pop", "pop"]], on="sample_id", how="left")
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_path, index=False)
    print(f"\nSaved {len(out):,} rows -> {out_path}")
    print(out.groupby("super_pop")["score_confidence"].mean())


if __name__ == "__main__":
    main()
