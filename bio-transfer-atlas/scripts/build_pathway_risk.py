"""
FAIRGEN-Open Stage 12: Pathway-level portability-risk aggregation (chr22)
==========================================================================
Maps chr22 variants -> genes (gene body +/-50kb) -> Reactome pathways, then
aggregates cross-ancestry concordance per (trait, pathway). Answers Q2:
which biological pathways drive cross-population score instability?

Inputs:
  data/raw/ensembl/chr22.gff3.gz                 gene coordinates (GRCh38)
  data/raw/reactome/Ensembl2Reactome_All_Levels.txt  gene -> pathway
  data/modeling/master_variant_table.parquet     labels + features per variant x trait

Outputs:
  data/annotations/variant_to_gene.parquet
  data/annotations/gene_to_pathway.parquet
  results/tables/pathway_risk_table.parquet
  results/tables/pathway_risk_top.csv
"""
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
E2R = ROOT / "data/raw/reactome/Ensembl2Reactome_All_Levels.txt"
ANN_DIR = ROOT / "data/annotations"
TAB_DIR = ROOT / "results/tables"
ANN_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

WINDOW = 50_000  # gene body +/- 50kb


def normalize_chr(value: str) -> str:
    s = str(value).strip()
    return s[3:] if s.lower().startswith("chr") else s


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22", help="Chromosome number/name (e.g. 22 or chr22)")
    return p.parse_args()


def resolve_gff_path(chrom: str) -> Path:
    p = ROOT / "data/raw/ensembl" / f"chr{chrom}.gff3.gz"
    if p.exists():
        return p
    if chrom == "22":
        legacy = ROOT / "data/raw/ensembl/chr22.gff3.gz"
        if legacy.exists():
            return legacy
    raise FileNotFoundError(f"Missing Ensembl GFF3 for chr{chrom}: {p}")


def resolve_master_path(chrom: str) -> Path:
    p = ROOT / "data/modeling" / f"master_variant_table.chr{chrom}.parquet"
    if p.exists():
        return p
    if chrom == "22":
        legacy = ROOT / "data/modeling/master_variant_table.parquet"
        if legacy.exists():
            return legacy
    raise FileNotFoundError(f"Missing master variant table for chr{chrom}: {p}")


def maybe_write_legacy(chrom: str, src: Path, legacy_name: str) -> None:
    if chrom == "22":
        if src.suffix == ".parquet":
            pd.read_parquet(src).to_parquet(src.parent / legacy_name, index=False)
        else:
            (src.parent / legacy_name).write_text(src.read_text(encoding="utf-8"), encoding="utf-8")


def main() -> None:
    args = parse_args()
    chrom = normalize_chr(args.chrom)
    gff = resolve_gff_path(chrom)
    master_p = resolve_master_path(chrom)

    print(f"Parsing chr{chrom} GFF3 ...")
    genes = []
    with gzip.open(gff, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            start, end, attrs = int(p[3]), int(p[4]), p[8]
            m_id = re.search(r"gene_id=([^;]+)", attrs)
            m_nm = re.search(r"Name=([^;]+)", attrs)
            m_bt = re.search(r"biotype=([^;]+)", attrs)
            if not m_id:
                continue
            genes.append(
                {
                    "ensg": m_id.group(1),
                    "gene_name": m_nm.group(1) if m_nm else "",
                    "start": start,
                    "end": end,
                    "biotype": m_bt.group(1) if m_bt else "",
                }
            )
    gdf = pd.DataFrame(genes)
    print(f"  {len(gdf):,} genes on chr{chrom}")

    master = pd.read_parquet(master_p)
    uv = master["variant_id"].drop_duplicates().to_frame()
    uv["pos"] = pd.to_numeric(uv["variant_id"].str.split(":").str[1], errors="coerce")
    uv = uv.dropna(subset=["pos"]).sort_values("pos").reset_index(drop=True)
    uv["pos"] = uv["pos"].astype(int)
    pos = uv["pos"].to_numpy()
    print(f"  {len(uv):,} unique variants")

    print("Mapping variants to genes (+/-50kb) ...")
    v2g_rows = []
    for g in gdf.itertuples(index=False):
        lo = np.searchsorted(pos, g.start - WINDOW, side="left")
        hi = np.searchsorted(pos, g.end + WINDOW, side="right")
        if hi > lo:
            for vid in uv["variant_id"].iloc[lo:hi]:
                v2g_rows.append((vid, g.ensg, g.gene_name))
    v2g = pd.DataFrame(v2g_rows, columns=["variant_id", "ensg", "gene_name"])
    v2g_path = ANN_DIR / f"variant_to_gene.chr{chrom}.parquet"
    v2g.to_parquet(v2g_path, index=False)
    maybe_write_legacy(chrom, v2g_path, "variant_to_gene.parquet")
    print(
        f"  {len(v2g):,} variant-gene links; {v2g['variant_id'].nunique():,} variants mapped to >=1 gene"
    )

    print("Loading Ensembl2Reactome ...")
    r = pd.read_csv(
        E2R,
        sep="\t",
        header=None,
        dtype=str,
        names=["source", "rid", "url", "pathway", "evidence", "species"],
    )
    r = r[(r["species"] == "Homo sapiens") & r["source"].str.startswith("ENSG")].copy()
    r = r[["source", "rid", "pathway"]].rename(columns={"source": "ensg"}).drop_duplicates()
    g2p_path = ANN_DIR / f"gene_to_pathway.chr{chrom}.parquet"
    r.to_parquet(g2p_path, index=False)
    maybe_write_legacy(chrom, g2p_path, "gene_to_pathway.parquet")
    print(f"  {len(r):,} gene-pathway links ({r['rid'].nunique():,} pathways)")

    vp = v2g.merge(r, on="ensg", how="inner", validate="many_to_many")[["variant_id", "rid", "pathway"]].drop_duplicates()
    print(f"  {vp['variant_id'].nunique():,} variants map to >=1 Reactome pathway")

    cols = [
        "variant_id",
        "trait",
        "I2",
        "sign_concordance",
        "risk_class",
        "associated",
        "PBS_max",
        "FST_like",
    ]
    mt = master[cols].merge(vp, on="variant_id", how="inner", validate="many_to_many")

    def agg(group: pd.DataFrame) -> pd.Series:
        a = group[group["associated"]]
        return pd.Series(
            {
                "n_variants": group["variant_id"].nunique(),
                "n_assoc": a["variant_id"].nunique(),
                "mean_I2": group["I2"].mean(),
                "mean_I2_assoc": a["I2"].mean() if len(a) else np.nan,
                "sign_discord_frac": (group["sign_concordance"] < 1.0).mean(),
                "high_risk_frac": (group["risk_class"] == 2).mean(),
                "mean_PBS": group["PBS_max"].mean(),
                "mean_FST": group["FST_like"].mean(),
            }
        )

    pr = mt.groupby(["trait", "rid", "pathway"]).apply(agg, include_groups=False).reset_index()
    pr = pr[pr["n_variants"] >= 5].copy()
    pr_path = TAB_DIR / f"pathway_risk_table.chr{chrom}.parquet"
    pr.to_parquet(pr_path, index=False)
    maybe_write_legacy(chrom, pr_path, "pathway_risk_table.parquet")
    print(f"\nSaved {pr_path.name} ({len(pr):,} trait-pathway rows)")

    lines = []

    def log(s: str) -> None:
        print(s)
        lines.append(s)

    log("=" * 70)
    log(f"TOP CROSS-ANCESTRY-UNSTABLE PATHWAYS per trait (chr{chrom})")
    log("=" * 70)
    top_all = []
    for trait in ["CAD", "T2D", "BMI", "LDL"]:
        sub = pr[(pr["trait"] == trait) & (pr["n_assoc"] >= 3)].copy()
        sub = sub.sort_values("mean_I2_assoc", ascending=False).head(8)
        if len(sub):
            top_all.append(sub)
        log(f"\n[{trait}]")
        if len(sub) == 0:
            log("  (no pathway with >=3 associated variants)")
            continue
        for x in sub.itertuples(index=False):
            log(
                f"  I2={x.mean_I2_assoc:.3f}  n_assoc={int(x.n_assoc):3d}  "
                f"PBS={x.mean_PBS:.3f}  {x.pathway[:60]}"
            )

    top_path = TAB_DIR / f"pathway_risk_top.chr{chrom}.csv"
    if top_all:
        pd.concat(top_all).to_csv(top_path, index=False, float_format="%.4f")
    else:
        pd.DataFrame(
            columns=["trait", "rid", "pathway", "n_variants", "n_assoc", "mean_I2_assoc", "mean_PBS"]
        ).to_csv(top_path, index=False)
    summary_path = TAB_DIR / f"pathway_risk_summary.chr{chrom}.txt"
    summary_path.write_text("\n".join(lines), encoding="utf-8")
    if chrom == "22":
        maybe_write_legacy(chrom, top_path, "pathway_risk_top.csv")
        maybe_write_legacy(chrom, summary_path, "pathway_risk_summary.txt")
    print(f"\nSaved {top_path.name} / {summary_path.name}")

if __name__ == "__main__":
    main()
