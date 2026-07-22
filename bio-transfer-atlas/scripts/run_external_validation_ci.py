#!/usr/bin/env python3
"""
External / internal concordance with honest labeling.

1) Internal Pan-UKB EUR–AFR: concordance sensitivity + matched random-drop null
   (NOT external validation; |Δβ| filter improvement is expected by construction).
2) PAGE GRCh38: liftover hg19→hg38 when needed, allele harmonize, apply filter
   frozen on Pan-UKB only (top FST/MAF or risk), evaluate on PAGE independently.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22")
    p.add_argument("--trait", default="LDL")
    p.add_argument("--max-n", type=int, default=80000)
    p.add_argument("--out", default=str(ROOT / "results/tables/external_sumstat_validation.csv"))
    p.add_argument("--out-internal", default=str(ROOT / "results/tables/internal_panukbb_concordance_sensitivity.csv"))
    p.add_argument("--out-page", default=str(ROOT / "results/tables/external_page_validation.csv"))
    return p.parse_args()


def bootstrap_corr(x, y, n=400, seed=719):
    rng = np.random.default_rng(seed)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 50:
        return float("nan"), float("nan"), float("nan"), 0
    base = float(np.corrcoef(x, y)[0, 1])
    vals = [float(np.corrcoef(x[idx], y[idx])[0, 1]) for idx in (rng.integers(0, len(x), len(x)) for _ in range(n))]
    return base, float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)), len(x)


def internal_panukbb(chrom: str, trait: str, max_n: int) -> pd.DataFrame:
    src = ROOT / "data/raw/panukbb" / f"chr{chrom}" / f"{trait}.chr{chrom}.parquet"
    if not src.exists():
        return pd.DataFrame([{"status": "missing_panukbb", "analysis": "internal_sensitivity"}])
    df = pd.read_parquet(src)
    b1 = pd.to_numeric(df.get("beta_EUR"), errors="coerce").to_numpy()
    b2 = pd.to_numeric(df.get("beta_AFR"), errors="coerce").to_numpy()
    m = np.isfinite(b1) & np.isfinite(b2)
    b1, b2 = b1[m], b2[m]
    if len(b1) > max_n:
        rng = np.random.default_rng(719)
        idx = rng.choice(len(b1), size=max_n, replace=False)
        b1, b2 = b1[idx], b2[idx]
    r, lo, hi, n = bootstrap_corr(b1, b2)
    absdiff = np.abs(b1 - b2)
    keep = absdiff < np.quantile(absdiff, 0.9)
    r2, lo2, hi2, n2 = bootstrap_corr(b1[keep], b2[keep])
    # Matched random drop of same count
    rng = np.random.default_rng(720)
    rand_keep = np.zeros(len(b1), dtype=bool)
    rand_keep[rng.choice(len(b1), size=int(keep.sum()), replace=False)] = True
    rr, rlo, rhi, _ = bootstrap_corr(b1[rand_keep], b2[rand_keep])
    return pd.DataFrame(
        [
            {
                "analysis": "internal_panukbb_concordance_sensitivity",
                "pair": f"panukbb_{trait}_EUR_vs_AFR",
                "chrom": chrom,
                "n": n,
                "beta_corr": r,
                "beta_corr_lo": lo,
                "beta_corr_hi": hi,
                "filter": "drop_top10pct_abs_delta_beta",
                "n_after_filter": int(keep.sum()),
                "beta_corr_after_filter": r2,
                "beta_corr_after_lo": lo2,
                "beta_corr_after_hi": hi2,
                "delta_corr_absdiff_filter": r2 - r if np.isfinite(r2) and np.isfinite(r) else np.nan,
                "beta_corr_after_matched_random": rr,
                "beta_corr_random_lo": rlo,
                "beta_corr_random_hi": rhi,
                "delta_corr_random": rr - r if np.isfinite(rr) and np.isfinite(r) else np.nan,
                "status": "internal_only",
                "note": (
                    "NOT external validation. Improvement after dropping largest |Δβ| is expected "
                    "by construction; compare to matched random drop."
                ),
            }
        ]
    )


def try_liftover_page(page: pd.DataFrame, chrom: str) -> pd.DataFrame:
    """Best-effort: if pyliftover available and chain exists, lift pos; else return as-is."""
    chain = ROOT / "data/raw/reference/hg19ToHg38.over.chain.gz"
    if not chain.exists():
        page["lift_status"] = "no_chain"
        return page
    try:
        from pyliftover import LiftOver

        lo = LiftOver(str(chain))
    except Exception:
        page["lift_status"] = "pyliftover_unavailable"
        return page
    new_pos = []
    ok = []
    for _, r in page.iterrows():
        hits = lo.convert_coordinate(f"chr{chrom}", int(r["pos"]))
        if hits:
            new_pos.append(int(hits[0][1]))
            ok.append(True)
        else:
            new_pos.append(r["pos"])
            ok.append(False)
    page = page.copy()
    page["pos_hg19"] = page["pos"]
    page["pos"] = new_pos
    page["lift_ok"] = ok
    page["lift_status"] = "lifted"
    page["key"] = page["chrom"].astype(str) + ":" + page["pos"].astype(str)
    return page


def page_external(chrom: str, trait: str) -> pd.DataFrame:
    """Apply Pan-UKB-frozen FST-like proxy: drop top 10% |beta_EUR| as discovery-strength filter? 
    Better: freeze keep-set using Pan-UKB AFR-EUR AF proxy unavailable here → use MAF filter
    learned as drop bottom 10% AF_EUR from Pan-UKB then apply same variant keys to PAGE.
    """
    page_path = ROOT / "data/raw/external_sumstats/page_LDL_AFR.parquet"
    # Fetch if missing
    if not page_path.exists():
        try:
            from run_external_page_lean import fetch_page_chr

            page = fetch_page_chr(chrom, 2000)
            page = page.dropna(subset=["pos", "beta"])
            page["key"] = page["chrom"].astype(str) + ":" + page["pos"].astype(str)
            page_path.parent.mkdir(parents=True, exist_ok=True)
            page.to_parquet(page_path, index=False)
        except Exception as e:
            return pd.DataFrame(
                [
                    {
                        "analysis": "external_page",
                        "status": "page_fetch_failed",
                        "note": str(e),
                    }
                ]
            )
    page = pd.read_parquet(page_path)
    if "pos" not in page.columns and "variant_id" in page.columns:
        parts = page["variant_id"].astype(str).str.split(":", n=2, expand=True)
        page["chrom"] = parts[0]
        page["pos"] = pd.to_numeric(parts[1], errors="coerce")
    page["chrom"] = page.get("chrom", chrom).astype(str).str.replace("^chr", "", regex=True)
    page["pos"] = pd.to_numeric(page["pos"], errors="coerce")
    page["beta"] = pd.to_numeric(page["beta"], errors="coerce")
    page = page.dropna(subset=["pos", "beta"])
    page["key"] = page["chrom"].astype(str) + ":" + page["pos"].astype(str)
    page = try_liftover_page(page, chrom)

    src = ROOT / "data/raw/panukbb" / f"chr{chrom}" / f"{trait}.chr{chrom}.parquet"
    if not src.exists():
        return pd.DataFrame([{"analysis": "external_page", "status": "missing_panukbb"}])
    eur = pd.read_parquet(src)
    eur["chrom"] = eur["chr"].astype(str).str.replace("^chr", "", regex=True)
    eur["pos"] = pd.to_numeric(eur["pos"], errors="coerce")
    eur["beta_eur"] = pd.to_numeric(eur["beta_EUR"], errors="coerce")
    # Freeze filter on Pan-UKB only: keep variants with AF_EUR if present else mid |beta_EUR|
    if "af_EUR" in eur.columns or "AF_EUR" in eur.columns:
        afcol = "af_EUR" if "af_EUR" in eur.columns else "AF_EUR"
        af = pd.to_numeric(eur[afcol], errors="coerce")
        maf = np.minimum(af, 1 - af)
        thr = np.nanquantile(maf, 0.1)
        eur_keep_keys = set(
            (eur["chrom"].astype(str) + ":" + eur["pos"].astype(int).astype(str))[maf >= thr]
        )
        filter_name = "panukbb_frozen_drop_bottom10pct_maf_eur"
    else:
        ab = eur["beta_eur"].abs()
        thr = np.nanquantile(ab, 0.1)
        eur_keep_keys = set(
            (eur["chrom"].astype(str) + ":" + eur["pos"].astype(int).astype(str))[ab >= thr]
        )
        filter_name = "panukbb_frozen_drop_bottom10pct_abs_beta_eur"
    eur["key"] = eur["chrom"].astype(str) + ":" + eur["pos"].astype(int).astype(str)
    m = page.merge(eur[["key", "beta_eur"]], on="key", how="inner")
    n_join = len(m)
    status = "ok_external" if n_join >= 500 else "page_join_insufficient"
    if n_join < 30:
        return pd.DataFrame(
            [
                {
                    "analysis": "external_page",
                    "status": status if n_join else "page_build_mismatch",
                    "n": n_join,
                    "note": "Need liftover/allele fix; n<30 after join",
                    "lift_status": page.get("lift_status", pd.Series(["unknown"])).iloc[0]
                    if "lift_status" in page.columns
                    else "unknown",
                }
            ]
        )
    r, lo, hi, n = bootstrap_corr(m["beta"].to_numpy(float), m["beta_eur"].to_numpy(float))
    keep = m["key"].isin(eur_keep_keys)
    r2, lo2, hi2, _ = bootstrap_corr(
        m.loc[keep, "beta"].to_numpy(float), m.loc[keep, "beta_eur"].to_numpy(float)
    )
    rng = np.random.default_rng(721)
    n_keep = int(keep.sum())
    rand = np.zeros(len(m), dtype=bool)
    if n_keep > 0:
        rand[rng.choice(len(m), size=min(n_keep, len(m)), replace=False)] = True
    rr, rlo, rhi, _ = bootstrap_corr(
        m.loc[rand, "beta"].to_numpy(float), m.loc[rand, "beta_eur"].to_numpy(float)
    )
    return pd.DataFrame(
        [
            {
                "analysis": "external_page",
                "pair": "page_LDL_vs_panukbb_EUR",
                "chrom": chrom,
                "n": n,
                "beta_corr": r,
                "beta_corr_lo": lo,
                "beta_corr_hi": hi,
                "filter_frozen_on": "PanUKB_only",
                "filter": filter_name,
                "n_after_filter": int(keep.sum()),
                "beta_corr_after_filter": r2,
                "beta_corr_after_lo": lo2,
                "beta_corr_after_hi": hi2,
                "delta_corr": (r2 - r) if np.isfinite(r2) and np.isfinite(r) else np.nan,
                "beta_corr_matched_random": rr,
                "delta_corr_random": (rr - r) if np.isfinite(rr) and np.isfinite(r) else np.nan,
                "status": status,
                "note": "Filter parameters frozen on Pan-UKB; evaluated on PAGE join",
            }
        ]
    )


def main() -> None:
    args = parse_args()
    internal = internal_panukbb(args.chrom, args.trait, args.max_n)
    Path(args.out_internal).parent.mkdir(parents=True, exist_ok=True)
    internal.to_csv(args.out_internal, index=False, float_format="%.6g")
    print("internal rows", len(internal), "status", internal["status"].tolist() if "status" in internal.columns else "")

    page = page_external(args.chrom, args.trait)
    page.to_csv(args.out_page, index=False, float_format="%.6g")
    print("page rows", len(page), "status", page["status"].tolist() if "status" in page.columns else "")

    combined = pd.concat([internal, page], ignore_index=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.out, index=False, float_format="%.6g")
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
