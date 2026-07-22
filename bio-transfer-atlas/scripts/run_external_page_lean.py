"""
Lean PAGE external validation via GWAS Catalog summary-statistics API (chr-filtered).

No FTP dumps. Pulls PAGE LDL associations for one chromosome and correlates
betas vs local Pan-UKB EUR LDL on chrom:pos.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/raw/external_sumstats"
GCST = "GCST008037"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22")
    p.add_argument("--size", type=int, default=1000)
    p.add_argument(
        "--out-table",
        default=str(ROOT / "results/tables/external_sumstat_validation.csv"),
    )
    return p.parse_args()


def bootstrap_corr(x: np.ndarray, y: np.ndarray, n: int = 200, seed: int = 719):
    rng = np.random.default_rng(seed)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 30:
        return float("nan"), float("nan"), float("nan")
    base = float(np.corrcoef(x, y)[0, 1])
    vals = [float(np.corrcoef(x[idx], y[idx])[0, 1]) for idx in (rng.integers(0, len(x), len(x)) for _ in range(n))]
    return base, float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def fetch_page_chr(chrom: str, size: int) -> pd.DataFrame:
    # Paginate summary-statistics API
    rows = []
    start = 0
    while len(rows) < size:
        url = (
            f"https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chrom}/associations"
            f"?study_accession={GCST}&size={min(500, size - len(rows))}&start={start}"
        )
        print(f"GET {url}", flush=True)
        r = requests.get(url, timeout=90)
        r.raise_for_status()
        data = r.json()
        assoc = data.get("_embedded", {}).get("associations", {})
        if isinstance(assoc, dict):
            batch = list(assoc.values())
        else:
            batch = assoc or []
        if not batch:
            break
        for a in batch:
            rows.append(
                {
                    "chrom": str(a.get("chromosome")),
                    "pos": int(a["base_pair_location"]) if a.get("base_pair_location") is not None else None,
                    "beta": a.get("beta"),
                    "se": None,
                    "pval": a.get("p_value"),
                    "effect_allele": a.get("effect_allele"),
                    "other_allele": a.get("other_allele"),
                    "variant_id_rs": a.get("variant_id"),
                }
            )
        start += len(batch)
        if len(batch) < 100:
            break
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    page = fetch_page_chr(args.chrom, args.size)
    page = page.dropna(subset=["pos", "beta"])
    page["pos"] = page["pos"].astype(int)
    page["beta"] = pd.to_numeric(page["beta"], errors="coerce")
    page = page.dropna(subset=["beta"])
    page["key"] = page["chrom"].astype(str) + ":" + page["pos"].astype(str)
    page_out = OUT / "page_LDL_AFR.parquet"
    page.rename(columns={"key": "variant_id"})[["variant_id", "beta", "pval", "chrom", "pos"]].to_parquet(
        page_out, index=False
    )
    print(f"Wrote {page_out} ({len(page):,})")

    src = ROOT / "data/raw/panukbb" / f"chr{args.chrom}" / f"LDL.chr{args.chrom}.parquet"
    if not src.exists():
        raise SystemExit(f"Missing {src}")
    import pyarrow.parquet as pq

    names = set(pq.read_schema(src).names)
    cols = [c for c in ["chr", "pos", "ref", "alt", "beta_EUR", "se_EUR", "pval_EUR", "pval_meta"] if c in names]
    eur = pd.read_parquet(src, columns=cols)
    eur["chrom"] = eur["chr"].astype(str).str.replace("^chr", "", regex=True)
    eur["pos"] = pd.to_numeric(eur["pos"], errors="coerce").astype("Int64")
    eur["beta"] = pd.to_numeric(eur["beta_EUR"], errors="coerce")
    eur["se"] = pd.to_numeric(eur["se_EUR"], errors="coerce") if "se_EUR" in eur.columns else np.nan
    pcol = "pval_EUR" if "pval_EUR" in eur.columns else ("pval_meta" if "pval_meta" in eur.columns else None)
    eur["pval"] = pd.to_numeric(eur[pcol], errors="coerce") if pcol else np.nan
    eur = eur.dropna(subset=["pos", "beta"])
    eur["key"] = eur["chrom"].astype(str) + ":" + eur["pos"].astype(str)
    eur_out = OUT / "panukbb_LDL_EUR.parquet"
    eur.assign(variant_id=eur["key"])[["variant_id", "beta", "se", "pval", "chrom", "pos"]].to_parquet(
        eur_out, index=False
    )

    m = page.merge(eur[["key", "beta"]], on="key", suffixes=("_page", "_eur"))
    print(f"Joined {len(m):,} PAGE x PanUKB-EUR", flush=True)
    r, lo, hi = bootstrap_corr(m["beta_page"].to_numpy(float), m["beta_eur"].to_numpy(float))
    thr = m["beta_page"].abs().quantile(0.9)
    m2 = m[m["beta_page"].abs() < thr]
    r2, lo2, hi2 = bootstrap_corr(m2["beta_page"].to_numpy(float), m2["beta_eur"].to_numpy(float))
    out = pd.DataFrame(
        [
            {
                "pair": "page_LDL_vs_panukbb_EUR",
                "anc": "PAGE_multi",
                "chrom": args.chrom,
                "n": len(m),
                "beta_corr": r,
                "beta_corr_lo": lo,
                "beta_corr_hi": hi,
                "beta_corr_after_top10pct_absbeta_filter": r2,
                "beta_corr_after_lo": lo2,
                "beta_corr_after_hi": hi2,
                "delta_corr": (r2 - r) if np.isfinite(r2) and np.isfinite(r) else np.nan,
                "status": "ok_lean_api",
                "note": "PAGE summary-statistics API chr-filtered; no FTP",
            }
        ]
    )
    Path(args.out_table).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_table, index=False, float_format="%.6g")
    print(out.to_string(index=False))
    print(f"Saved {args.out_table}")


if __name__ == "__main__":
    main()
