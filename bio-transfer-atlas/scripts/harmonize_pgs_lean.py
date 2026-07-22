#!/usr/bin/env python3
"""
Lean harmonize for small expansion PGS (WBC/RA/IBD): stream pvars for needed positions only.
Does not rebuild full-chromosome lookup caches.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCORES_RAW = ROOT / "data/raw/pgs_catalog/scores"
PGS_OUT = ROOT / "data/processed/pgs_grch38"
INTERIM = ROOT / "data/interim/1000g_grch38"
AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

DEFAULT_IDS = ["PGS000191", "PGS004133", "PGS001288"]


def _find(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def read_pgs(pgs_id: str) -> pd.DataFrame | None:
    src = SCORES_RAW / pgs_id / f"{pgs_id}_hmPOS_GRCh38.txt.gz"
    if not src.exists():
        print(f"missing {src}")
        return None
    df = pd.read_csv(src, sep="\t", compression="gzip", comment="#", dtype=str, low_memory=False)
    cols = df.columns.tolist()
    rename = {}
    for src_c, dst in [
        (_find(cols, ["hm_chr", "chr_name"]), "chr"),
        (_find(cols, ["hm_pos", "chr_position"]), "pos"),
        (_find(cols, ["effect_allele"]), "effect_allele"),
        (_find(cols, ["other_allele"]), "other_allele"),
        (_find(cols, ["effect_weight"]), "effect_weight"),
    ]:
        if src_c:
            rename[src_c] = dst
    df = df[list(rename)].rename(columns=rename)
    df["chr"] = df["chr"].astype(str).str.replace("^chr", "", regex=True)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df["effect_weight"] = pd.to_numeric(df["effect_weight"], errors="coerce")
    df = df.dropna(subset=["chr", "pos", "effect_allele", "effect_weight"])
    return df


def stream_match(chrom: str, need_pos: set[int]) -> pd.DataFrame:
    pvar = INTERIM / f"chr{chrom}.score.pvar"
    if not pvar.exists() or not need_pos:
        return pd.DataFrame(columns=["chr", "pos", "ref", "alt"])
    rows = []
    with pvar.open("r", encoding="utf-8", errors="replace") as f:
        header = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#") and header is None:
                header = line.lstrip("#").strip().split("\t")
                continue
            parts = line.rstrip("\n").split("\t")
            # PLINK2 pvar: CHROM POS ID REF ALT ...
            try:
                pos = int(float(parts[1]))
            except Exception:
                continue
            if pos not in need_pos:
                continue
            chrom_v = parts[0].replace("chr", "")
            ref, alt = parts[3], parts[4].split(",")[0]
            rows.append({"chr": chrom_v, "pos": pos, "ref": ref, "alt": alt})
    return pd.DataFrame(rows)


def match_alleles(pgs_sub: pd.DataFrame, pvar_sub: pd.DataFrame) -> pd.DataFrame:
    if pvar_sub.empty or pgs_sub.empty:
        return pd.DataFrame(columns=["variant_id", "effect_allele", "effect_weight"])
    pvar_idx = pvar_sub.drop_duplicates(["chr", "pos"]).set_index(["chr", "pos"])
    m = pgs_sub.set_index(["chr", "pos"]).join(pvar_idx[["ref", "alt"]], how="inner").reset_index()
    if m.empty:
        return pd.DataFrame(columns=["variant_id", "effect_allele", "effect_weight"])
    keep = []
    for _, r in m.iterrows():
        ea = str(r["effect_allele"]).upper()
        ref, alt = str(r["ref"]).upper(), str(r["alt"]).upper()
        if (ea, ref) in AMBIGUOUS or (ea, alt) in AMBIGUOUS:
            # keep if matches ref or alt
            pass
        if ea not in (ref, alt):
            continue
        vid = f"{r['chr']}:{int(r['pos'])}:{ref}:{alt}"
        keep.append({"variant_id": vid, "effect_allele": ea, "effect_weight": r["effect_weight"]})
    return pd.DataFrame(keep)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pgs-ids", default=",".join(DEFAULT_IDS))
    args = ap.parse_args()
    ids = [x.strip() for x in args.pgs_ids.split(",") if x.strip()]

    for pgs_id in ids:
        out = PGS_OUT / pgs_id / f"{pgs_id}.harmonized.tsv"
        if out.exists() and out.stat().st_size > 50:
            print(f"{pgs_id}: already harmonized ({out.stat().st_size} B)")
            continue
        pgs = read_pgs(pgs_id)
        if pgs is None:
            continue
        print(f"{pgs_id}: {len(pgs)} variants", flush=True)
        frames = []
        for chrom, sub in pgs.groupby(pgs["chr"].astype(str)):
            need = set(int(x) for x in sub["pos"].tolist())
            print(f"  chr{chrom}: scan {len(need)} positions", flush=True)
            pvar_sub = stream_match(str(chrom), need)
            hits = match_alleles(sub, pvar_sub)
            if not hits.empty:
                frames.append(hits)
                print(f"    matched {len(hits)}", flush=True)
        result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
            columns=["variant_id", "effect_allele", "effect_weight"]
        )
        out.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(out, sep="\t", index=False)
        print(f"  wrote {out} ({len(result)} rows)", flush=True)
    print("LEAN_HARMONIZE_DONE")


if __name__ == "__main__":
    main()
