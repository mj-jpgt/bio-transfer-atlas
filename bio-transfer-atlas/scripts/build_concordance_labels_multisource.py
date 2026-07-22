"""
Build multi-source cross-ancestry concordance labels.

Sources:
- Pan-UKBB: AFR/AMR/CSA/EAS/EUR/MID (GRCh37)
- BBJ: EAS (GRCh37)
- FinnGen: EUR-FIN (GRCh38)
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from pyliftover import LiftOver
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
CHAIN = ROOT / "data/raw/reference/hg38ToHg19.over.chain.gz"
OUT_DIR = ROOT / "data/labels"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TRAITS = ["T2D", "CAD", "BMI", "LDL"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22", help="Chromosome number/name (e.g. 22 or chr22)")
    return p.parse_args()


def normalize_chr(value: str) -> str:
    s = str(value).strip()
    return s[3:] if s.lower().startswith("chr") else s


def resolve_af_path(chrom: str) -> Path:
    af = ROOT / "data/features/af" / f"af_divergence_features.chr{chrom}.parquet"
    if af.exists():
        return af
    if chrom == "22":
        legacy = ROOT / "data/features/af/af_divergence_features.parquet"
        if legacy.exists():
            return legacy
    raise FileNotFoundError(
        f"Missing AF features for chr{chrom}. Run compute_af_features.py --chrom {chrom} first."
    )


def maybe_write_legacy(chrom: str, src: Path, legacy_name: str) -> None:
    if chrom == "22":
        pd.read_parquet(src).to_parquet(OUT_DIR / legacy_name, index=False)


def allele_key(ref: pd.Series, alt: pd.Series) -> pd.Series:
    a = ref.str.upper()
    b = alt.str.upper()
    return np.where(a <= b, a + "|" + b, b + "|" + a)


def build_target_maps(chrom: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    feat = pd.read_parquet(resolve_af_path(chrom), columns=["variant_id"])
    p = feat["variant_id"].str.split(":", expand=True)
    feat["chr"] = p[0].astype(str)
    feat["pos38"] = pd.to_numeric(p[1], errors="coerce").astype("Int64")
    feat["ref"] = p[2].str.upper()
    feat["alt"] = p[3].str.upper()
    feat = feat.dropna(subset=["pos38"]).copy()
    feat = feat[feat["chr"].map(normalize_chr) == chrom].copy()
    feat["pos38"] = feat["pos38"].astype(int)
    feat["allele_key"] = allele_key(feat["ref"], feat["alt"])
    map38 = feat[["variant_id", "pos38", "allele_key"]].drop_duplicates()

    cache = OUT_DIR / f"_chr{chrom}_pos37_cache.parquet"
    if cache.exists():
        c = pd.read_parquet(cache)
        pos37 = dict(zip(c["pos38"], c["pos37"]))
    else:
        lo = LiftOver(str(CHAIN))
        pos37 = {}
        for pos in map38["pos38"].unique():
            conv = lo.convert_coordinate(f"chr{chrom}", int(pos))
            if conv:
                mapped_chrom, new_pos, *_ = conv[0]
                if normalize_chr(mapped_chrom) == chrom:
                    pos37[int(pos)] = int(new_pos)
        pd.DataFrame({"pos38": list(pos37), "pos37": list(pos37.values())}).to_parquet(cache, index=False)
    map37 = map38.copy()
    map37["pos37"] = map37["pos38"].map(pos37)
    map37 = map37.dropna(subset=["pos37"]).copy()
    map37["pos37"] = map37["pos37"].astype(int)
    map37 = map37[["variant_id", "pos37", "allele_key"]].drop_duplicates()
    return map37, map38


def from_pan(map37: pd.DataFrame, chrom: str, pan_dir: Path, trait_filter: str | None = None) -> pd.DataFrame:
    rows: List[pd.DataFrame] = []
    for trait in ([trait_filter] if trait_filter else TRAITS):
        fp = pan_dir / f"{trait}.chr{chrom}.parquet"
        if not fp.exists():
            continue
        d = pd.read_parquet(fp)
        d["pos37"] = pd.to_numeric(d["pos"], errors="coerce")
        d = d.dropna(subset=["pos37"]).copy()
        d["pos37"] = d["pos37"].astype(int)
        d["allele_key"] = allele_key(d["ref"].astype(str), d["alt"].astype(str))
        m = d.merge(map37, on=["pos37", "allele_key"], how="inner", validate="many_to_one")
        anc_cols = [c.replace("beta_", "") for c in m.columns if c.startswith("beta_")]
        for anc in anc_cols:
            b = f"beta_{anc}"
            s = f"se_{anc}"
            if s not in m.columns:
                continue
            sub = m[["variant_id", b, s]].copy()
            sub = sub.rename(columns={b: "beta", s: "se"})
            sub["trait"] = trait
            sub["source"] = f"PanUKBB_{anc}"
            rows.append(sub)
    if not rows:
        return pd.DataFrame(columns=["variant_id", "trait", "source", "beta", "se"])
    out = pd.concat(rows, ignore_index=True)
    out["beta"] = pd.to_numeric(out["beta"], errors="coerce")
    out["se"] = pd.to_numeric(out["se"], errors="coerce")
    out = out.dropna(subset=["beta", "se"])
    out = out[out["se"] > 0]
    return out


def from_bbj(map37: pd.DataFrame, chrom: str, bbj_dir: Path, trait_filter: str | None = None) -> pd.DataFrame:
    rows = []
    for trait in ([trait_filter] if trait_filter else TRAITS):
        fp = bbj_dir / f"{trait}.chr{chrom}.parquet"
        if not fp.exists():
            continue
        d = pd.read_parquet(fp)
        d["pos37"] = pd.to_numeric(d["pos"], errors="coerce")
        d = d.dropna(subset=["pos37"]).copy()
        d["pos37"] = d["pos37"].astype(int)
        d["allele_key"] = allele_key(d["ref"].astype(str), d["alt"].astype(str))
        m = d.merge(map37, on=["pos37", "allele_key"], how="inner", validate="many_to_one")
        sub = m[["variant_id", "beta", "se"]].copy()
        sub["trait"] = trait
        sub["source"] = "BBJ_EAS"
        rows.append(sub)
    if not rows:
        return pd.DataFrame(columns=["variant_id", "trait", "source", "beta", "se"])
    out = pd.concat(rows, ignore_index=True)
    out["beta"] = pd.to_numeric(out["beta"], errors="coerce")
    out["se"] = pd.to_numeric(out["se"], errors="coerce")
    out = out.dropna(subset=["beta", "se"])
    out = out[out["se"] > 0]
    return out


def from_finngen(map38: pd.DataFrame, chrom: str, fin_dir: Path, trait_filter: str | None = None) -> pd.DataFrame:
    # FinnGen currently contributes CAD + T2D disease endpoints.
    rows = []
    finngen_traits = ["T2D", "CAD"]
    for trait in ([trait_filter] if (trait_filter and trait_filter in finngen_traits) else finngen_traits):
        fp = fin_dir / f"{trait}.chr{chrom}.parquet"
        if not fp.exists():
            continue
        d = pd.read_parquet(fp)
        d["pos38"] = pd.to_numeric(d["pos"], errors="coerce")
        d = d.dropna(subset=["pos38"]).copy()
        d["pos38"] = d["pos38"].astype(int)
        d["allele_key"] = allele_key(d["ref"].astype(str), d["alt"].astype(str))
        m = d.merge(map38, on=["pos38", "allele_key"], how="inner", validate="many_to_one")
        sub = m[["variant_id", "beta", "se"]].copy()
        sub["trait"] = trait
        sub["source"] = "FinnGen_EUR_FIN"
        rows.append(sub)
    if not rows:
        return pd.DataFrame(columns=["variant_id", "trait", "source", "beta", "se"])
    out = pd.concat(rows, ignore_index=True)
    out["beta"] = pd.to_numeric(out["beta"], errors="coerce")
    out["se"] = pd.to_numeric(out["se"], errors="coerce")
    out = out.dropna(subset=["beta", "se"])
    out = out[out["se"] > 0]
    return out


def concordance_metrics(g: pd.DataFrame) -> Dict:
    betas = g["beta"].to_numpy(float)
    ses = g["se"].to_numpy(float)
    srcs = g["source"].astype(str).tolist()
    n = len(betas)
    w = 1.0 / (ses ** 2)
    bmean = np.average(betas, weights=w)
    se_fe = np.sqrt(1.0 / np.sum(w))
    z_meta = bmean / se_fe if se_fe > 0 else np.nan
    Q = float(np.sum(w * (betas - bmean) ** 2))
    dfree = n - 1
    I2 = max(0.0, (Q - dfree) / Q) if Q > 0 else 0.0
    het_p = float(stats.chi2.sf(Q, dfree)) if dfree > 0 else np.nan
    signs = np.sign(betas)
    pairs = list(combinations(range(n), 2))
    sign_conc = float(np.mean([signs[i] == signs[j] for i, j in pairs])) if pairs else 1.0
    src_pair = {}
    for i, j in pairs:
        key = f"{srcs[i]}__{srcs[j]}"
        src_pair[key] = int(signs[i] == signs[j])
    risk = 0
    if I2 > 0.25 or sign_conc < 0.80:
        risk = 1
    if I2 > 0.50 or sign_conc < 0.60:
        risk = 2
    only_pan = all(s.startswith("PanUKBB_") for s in srcs)
    return {
        "n_sources": n,
        "sources": ",".join(srcs),
        "source_pair_sign_json": json.dumps(src_pair, separators=(",", ":")),
        "beta_fixed": float(bmean),
        "se_fixed": float(se_fe),
        "z_meta": float(z_meta) if pd.notna(z_meta) else np.nan,
        "associated": bool(abs(z_meta) >= 3.0) if pd.notna(z_meta) else False,
        "sign_concordance": sign_conc,
        "cochran_Q": Q,
        "I2": float(I2),
        "het_pval": het_p,
        "risk_class": risk,
        "source_group": "panukbb_only" if only_pan else "multisource",
    }


if __name__ == "__main__":
    args = parse_args()
    chrom = normalize_chr(args.chrom)
    pan_dir = ROOT / "data/raw/panukbb" / f"chr{chrom}"
    bbj_dir = ROOT / "data/raw/bbj" / f"chr{chrom}"
    fin_dir = ROOT / "data/raw/finngen" / f"chr{chrom}"

    print(f"Preparing target maps for chr{chrom} ...")
    map37, map38 = build_target_maps(chrom)
    print(f"target map37={len(map37):,} map38={len(map38):,}")

    # Process one trait at a time; write directly to parquet to cap peak memory
    import pyarrow as pa
    import pyarrow.parquet as pq
    out_path = OUT_DIR / f"gwas_concordance_labels_multisource.chr{chrom}.parquet"
    pq_writer: pq.ParquetWriter | None = None
    total_source_rows = 0
    for trait in TRAITS:
        pan = from_pan(map37, chrom, pan_dir, trait_filter=trait)
        bbj = from_bbj(map37, chrom, bbj_dir, trait_filter=trait)
        fin = from_finngen(map38, chrom, fin_dir, trait_filter=trait)
        trait_src = pd.concat([pan, bbj, fin], ignore_index=True)
        del pan, bbj, fin
        trait_src = trait_src.dropna(subset=["variant_id", "trait", "beta", "se"]).copy()
        total_source_rows += len(trait_src)
        print(f"  {trait}: {len(trait_src):,} source rows | " +
              trait_src["source"].value_counts().to_string().replace("\n", ", "), flush=True)

        trait_rows = []
        for variant_id, g in trait_src.groupby("variant_id", sort=False):
            if len(g) < 2:
                continue
            m = concordance_metrics(g)
            m["variant_id"] = variant_id
            m["trait"] = trait
            trait_rows.append(m)
        del trait_src

        if trait_rows:
            tf = pd.DataFrame(trait_rows)
            del trait_rows
            tf = tf.sort_values("variant_id").reset_index(drop=True)
            tf["associated_multisource"] = tf["associated"] & (tf["source_group"] == "multisource")
            for col in tf.select_dtypes(include="float64").columns:
                tf[col] = tf[col].astype("float32")
            # Write this trait directly to parquet (avoids large concat later)
            table = pa.Table.from_pandas(tf, preserve_index=False)
            del tf
            if pq_writer is None:
                pq_writer = pq.ParquetWriter(str(out_path), table.schema)
            pq_writer.write_table(table)
            del table
            print(f"  {trait}: wrote to parquet", flush=True)

    if pq_writer is not None:
        pq_writer.close()

    print(f"source rows total: {total_source_rows:,}")
    maybe_write_legacy(chrom, out_path, "gwas_concordance_labels_multisource.parquet")

    # Read back for audit (parquet is already compact; reading is manageable)
    labels = pd.read_parquet(out_path)
    print(f"Saved {out_path} with {len(labels):,} rows")

    lines = []

    def log(s: str):
        print(s)
        lines.append(s)

    log("=" * 64)
    log("MULTI-SOURCE LABEL AUDIT")
    log("=" * 64)
    log(f"Rows: {len(labels):,}")
    log(f"Unique variants: {labels['variant_id'].nunique():,}")
    log("Rows by trait:")
    log(labels["trait"].value_counts().to_string())
    log("\nRows by source_group:")
    log(labels["source_group"].value_counts().to_string())
    log("\nn_sources distribution:")
    log(labels["n_sources"].value_counts().sort_index().to_string())
    assoc = labels[labels["associated"]]
    assoc_ms = labels[labels["associated_multisource"]]
    log(f"\nAssociated rows: {len(assoc):,} ({len(assoc)/max(1,len(labels))*100:.1f}%)")
    if len(assoc):
        log(f"Associated median I2: {assoc['I2'].median():.4f}")
    log(f"Associated+multisource rows: {len(assoc_ms):,}")
    if len(assoc_ms):
        log(f"Associated+multisource median I2: {assoc_ms['I2'].median():.4f}")
    audit_path = OUT_DIR / f"label_audit_multisource.chr{chrom}.txt"
    audit_path.write_text("\n".join(lines), encoding="utf-8")
    if chrom == "22":
        (OUT_DIR / "label_audit_multisource.txt").write_text("\n".join(lines), encoding="utf-8")
