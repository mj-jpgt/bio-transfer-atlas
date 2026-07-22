#!/usr/bin/env python3
"""
PAGE GRCh38 external validation with allele QC ladder + frozen Pan-UKB filter.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/raw/external_sumstats"
GCST = "GCST008037"
TABLES = ROOT / "results/tables"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chroms", default="1,2,6,19,22")
    p.add_argument("--size-per-chrom", type=int, default=20000)
    p.add_argument("--trait", default="LDL")
    p.add_argument("--out", default=str(TABLES / "external_page_validation.csv"))
    p.add_argument(
        "--from-joined",
        default="",
        help="Reuse cached PAGE×Pan-UKB join parquet (skip API fetch)",
    )
    return p.parse_args()


def bootstrap_corr(x, y, n=400, seed=719):
    rng = np.random.default_rng(seed)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[mask], np.asarray(y)[mask]
    if len(x) < 50:
        return float("nan"), float("nan"), float("nan"), 0
    base = float(np.corrcoef(x, y)[0, 1])
    vals = [
        float(np.corrcoef(x[idx], y[idx])[0, 1])
        for idx in (rng.integers(0, len(x), len(x)) for _ in range(n))
    ]
    return base, float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)), len(x)


def fetch_page_chr(chrom: str, size: int) -> pd.DataFrame:
    rows = []
    start = 0
    while len(rows) < size:
        url = (
            f"https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chrom}/associations"
            f"?study_accession={GCST}&size={min(500, size - len(rows))}&start={start}"
        )
        print(f"GET {url}", flush=True)
        data = None
        for attempt in range(8):
            r = requests.get(url, timeout=120)
            if r.status_code == 429:
                wait = min(120, 5 * (2**attempt))
                print(f"  429; sleep {wait}s", flush=True)
                time.sleep(wait)
                continue
            r.raise_for_status()
            data = r.json()
            break
        if data is None:
            raise RuntimeError(f"PAGE API rate-limited on chr{chrom} start={start}")
        assoc = data.get("_embedded", {}).get("associations", {})
        batch = list(assoc.values()) if isinstance(assoc, dict) else (assoc or [])
        if not batch:
            break
        for a in batch:
            rows.append(
                {
                    "chrom": str(a.get("chromosome")),
                    "pos_hg19": int(a["base_pair_location"])
                    if a.get("base_pair_location") is not None
                    else None,
                    "beta": a.get("beta"),
                    "pval": a.get("p_value"),
                    "effect_allele": a.get("effect_allele"),
                    "other_allele": a.get("other_allele"),
                    "rsid": a.get("variant_id"),
                }
            )
        start += len(batch)
        time.sleep(0.35)
        if len(batch) < 100:
            break
    return pd.DataFrame(rows)


def liftover(df: pd.DataFrame) -> pd.DataFrame:
    chain = ROOT / "data/raw/reference/hg19ToHg38.over.chain.gz"
    if not chain.exists():
        df = df.copy()
        df["pos"] = df["pos_hg19"]
        df["lift_status"] = "no_chain"
        return df
    try:
        from pyliftover import LiftOver

        lo = LiftOver(str(chain))
    except Exception:
        df = df.copy()
        df["pos"] = df["pos_hg19"]
        df["lift_status"] = "pyliftover_unavailable"
        return df
    pos38, ok = [], []
    for _, r in df.iterrows():
        chrom = str(r["chrom"]).replace("chr", "")
        try:
            hits = lo.convert_coordinate(f"chr{chrom}", int(r["pos_hg19"]))
        except Exception:
            hits = None
        if hits:
            pos38.append(int(hits[0][1]))
            ok.append(True)
        else:
            pos38.append(-1)
            ok.append(False)
    out = df.copy()
    out["pos"] = pos38
    out["lift_ok"] = ok
    out["lift_status"] = "lifted"
    return out


def load_panukbb(chrom: str, trait: str) -> pd.DataFrame:
    src = ROOT / "data/raw/panukbb" / f"chr{chrom}" / f"{trait}.chr{chrom}.parquet"
    if not src.exists():
        return pd.DataFrame()
    eur = pd.read_parquet(src)
    eur["chrom"] = eur["chr"].astype(str).str.replace("^chr", "", regex=True)
    eur["pos"] = pd.to_numeric(eur["pos"], errors="coerce").astype("Int64")
    eur["beta_eur"] = pd.to_numeric(eur["beta_EUR"], errors="coerce")
    for c in ("ref", "alt", "allele1", "allele2", "REF", "ALT"):
        if c in eur.columns:
            eur[c] = eur[c].astype(str)
    eur = eur.dropna(subset=["pos", "beta_eur"])
    eur["key"] = eur["chrom"].astype(str) + ":" + eur["pos"].astype(int).astype(str)
    return eur


def freeze_keep_keys(eur: pd.DataFrame) -> tuple[set[str], str, list[str]]:
    """Learn keep-set on Pan-UKB only. Returns keys, filter name, training columns used."""
    train_cols = ["beta_eur"]
    if "af_EUR" in eur.columns or "AF_EUR" in eur.columns:
        afcol = "af_EUR" if "af_EUR" in eur.columns else "AF_EUR"
        train_cols.append(afcol)
        af = pd.to_numeric(eur[afcol], errors="coerce")
        maf = np.minimum(af, 1.0 - af)
        thr = np.nanquantile(maf, 0.1)
        keep = set(eur.loc[maf >= thr, "key"].astype(str))
        return keep, "panukbb_frozen_drop_bottom10pct_maf_eur", train_cols
    ab = eur["beta_eur"].abs()
    thr = np.nanquantile(ab, 0.1)
    keep = set(eur.loc[ab >= thr, "key"].astype(str))
    return keep, "panukbb_frozen_drop_bottom10pct_abs_beta_eur", train_cols


def _is_ambiguous(a: str, b: str) -> bool:
    a, b = a.upper(), b.upper()
    return {a, b} in ({"A", "T"}, {"C", "G"})


def allele_harmonize(m: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Match / swap+flip / complement / complement+swap; drop ambiguous & indels."""
    counts = {
        "position_joined": len(m),
        "allele_harmonized": 0,
        "ambiguous_removed": 0,
        "indel_removed": 0,
        "unresolved_removed": 0,
    }
    if "effect_allele" not in m.columns:
        m = m.copy()
        m["allele_status"] = "no_page_alleles"
        counts["allele_harmonized"] = len(m)
        return m, counts
    ref = next((c for c in ("ref", "allele1", "REF") if c in m.columns), None)
    alt = next((c for c in ("alt", "allele2", "ALT") if c in m.columns), None)
    if ref is None or alt is None:
        m = m.copy()
        m["allele_status"] = "no_panukbb_alleles"
        counts["allele_harmonized"] = len(m)
        return m, counts

    rows = []
    for _, r in m.iterrows():
        ea = str(r.get("effect_allele", "")).upper()
        oa = str(r.get("other_allele", "")).upper()
        rr = str(r[ref]).upper()
        aa = str(r[alt]).upper()
        beta = float(r["beta"]) if pd.notna(r["beta"]) else np.nan
        if len(ea) != 1 or len(oa) != 1 or len(rr) != 1 or len(aa) != 1:
            counts["indel_removed"] += 1
            continue
        if _is_ambiguous(ea, oa) or _is_ambiguous(rr, aa):
            counts["ambiguous_removed"] += 1
            continue
        status = None
        b = beta
        if ea == aa and oa == rr:
            status = "direct"
        elif ea == rr and oa == aa:
            status = "swap_flip"
            b = -beta
        elif COMPLEMENT.get(ea) == aa and COMPLEMENT.get(oa) == rr:
            status = "complement"
        elif COMPLEMENT.get(ea) == rr and COMPLEMENT.get(oa) == aa:
            status = "complement_swap_flip"
            b = -beta
        else:
            counts["unresolved_removed"] += 1
            continue
        out = r.to_dict()
        out["beta"] = b
        out["allele_status"] = status
        rows.append(out)
        counts["allele_harmonized"] += 1
    return pd.DataFrame(rows), counts


def richer_metrics(x: np.ndarray, y: np.ndarray) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 50:
        return {
            "pearson_r": np.nan,
            "spearman_r": np.nan,
            "sign_concordance": np.nan,
            "regression_slope": np.nan,
            "standardized_beta_rmse": np.nan,
        }
    pearson = float(np.corrcoef(x, y)[0, 1])
    from scipy import stats

    try:
        spearman = float(stats.spearmanr(x, y).correlation)
    except Exception:
        spearman = float("nan")
    sign_conc = float(np.mean(np.sign(x) == np.sign(y)))
    slope = float(np.polyfit(x, y, 1)[0])
    xz = (x - x.mean()) / (x.std() + 1e-12)
    yz = (y - y.mean()) / (y.std() + 1e-12)
    rmse = float(np.sqrt(np.mean((xz - yz) ** 2)))
    return {
        "pearson_r": pearson,
        "spearman_r": spearman,
        "sign_concordance": sign_conc,
        "regression_slope": slope,
        "standardized_beta_rmse": rmse,
    }


def main() -> None:
    args = parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    qc = {
        "input_page": 0,
        "lifted": 0,
        "position_joined": 0,
        "allele_harmonized": 0,
        "ambiguous_removed": 0,
        "indel_removed": 0,
        "unresolved_removed": 0,
        "final_variants": 0,
    }
    parts = []
    filter_name = "unset"
    train_cols: list[str] = []

    cache = Path(args.from_joined) if args.from_joined else OUT / "page_panukbb_joined_grch38.parquet"
    use_cache = bool(args.from_joined) or cache.exists()
    if use_cache:
        # Prefer cache to avoid GWAS Catalog rate limits; re-apply allele QC ladder
        print(f"Using cached join {cache}", flush=True)
        m0 = pd.read_parquet(cache)
        qc["input_page"] = len(m0)
        qc["lifted"] = int(m0["pos"].gt(0).sum()) if "pos" in m0.columns else len(m0)
        page_cols = {"beta", "effect_allele", "other_allele", "rsid", "pos_hg19"}
        train_cols = ["beta_eur"]
        assert page_cols.isdisjoint(set(train_cols))
        if "eur_keep" not in m0.columns or "filter_name" not in m0.columns:
            keep_all = set()
            for chrom in sorted(m0["chrom"].astype(str).unique()):
                eur = load_panukbb(str(chrom), args.trait)
                if eur.empty:
                    continue
                keep_keys, filter_name, train_cols = freeze_keep_keys(eur)
                keep_all |= keep_keys
            m0["eur_keep"] = m0["key"].astype(str).isin(keep_all)
            m0["filter_name"] = filter_name
        else:
            filter_name = str(m0["filter_name"].iloc[0])
        m, ac = allele_harmonize(m0)
        for k, v in ac.items():
            qc[k] = qc.get(k, 0) + v
        if not m.empty:
            parts = [m]
            print(f"cache allele-QC: n={len(m)} ambiguous={qc.get('ambiguous_removed', 0)}", flush=True)
    else:
        for chrom in chroms:
            page = fetch_page_chr(chrom, args.size_per_chrom)
            if page.empty:
                continue
            page = page.dropna(subset=["pos_hg19", "beta"])
            page["beta"] = pd.to_numeric(page["beta"], errors="coerce")
            page = page.dropna(subset=["beta"])
            qc["input_page"] += len(page)
            page_cols = set(page.columns)
            page = liftover(page)
            page = page[page["pos"] > 0].copy()
            qc["lifted"] += len(page)
            page["key"] = page["chrom"].astype(str) + ":" + page["pos"].astype(int).astype(str)
            eur = load_panukbb(chrom, args.trait)
            if eur.empty:
                continue
            keep_keys, filter_name, train_cols = freeze_keep_keys(eur)
            assert page_cols.isdisjoint(set(train_cols)), "PAGE columns used to train filter"
            cols = ["key", "beta_eur"] + [
                c for c in ("ref", "alt", "allele1", "allele2", "REF", "ALT") if c in eur.columns
            ]
            m = page.merge(eur[cols], on="key", how="inner")
            m, ac = allele_harmonize(m)
            for k in (
                "position_joined",
                "allele_harmonized",
                "ambiguous_removed",
                "indel_removed",
                "unresolved_removed",
            ):
                qc[k] = qc.get(k, 0) + ac.get(k, 0)
            if m.empty:
                continue
            m["chrom"] = chrom
            m["filter_name"] = filter_name
            m["eur_keep"] = m["key"].isin(keep_keys)
            parts.append(m)
            print(f"chr{chrom}: lifted={qc['lifted']:,} joined_harm={len(m):,}", flush=True)

    if not parts:
        out = pd.DataFrame(
            [
                {
                    "analysis": "external_page",
                    "status": "page_join_insufficient",
                    "n_variants": 0,
                    **{f"qc_{k}": v for k, v in qc.items()},
                }
            ]
        )
        out.to_csv(args.out, index=False)
        print(out.to_string(index=False))
        return

    m = pd.concat(parts, ignore_index=True)
    # drop duplicate keys
    m = m.drop_duplicates("key")
    qc["final_variants"] = len(m)
    m.to_parquet(OUT / "page_panukbb_joined_grch38.parquet", index=False)

    x = m["beta"].to_numpy(float)
    y = m["beta_eur"].to_numpy(float)
    metrics = richer_metrics(x, y)
    r, lo, hi, n = bootstrap_corr(x, y)
    keep = m["eur_keep"].to_numpy(bool)
    r2, lo2, hi2, n2 = bootstrap_corr(x[keep], y[keep])
    rng = np.random.default_rng(721)
    rand = np.zeros(len(m), dtype=bool)
    n_keep = int(keep.sum())
    if n_keep > 0:
        rand[rng.choice(len(m), size=min(n_keep, len(m)), replace=False)] = True
    rr, _, _, _ = bootstrap_corr(x[rand], y[rand])
    status = "ok_external" if n >= 500 else "page_join_insufficient"
    out = pd.DataFrame(
        [
            {
                "analysis": "external_page",
                "pair": "page_LDL_vs_panukbb_EUR",
                "chroms": ",".join(chroms),
                "n_variants": n,
                "n": n,  # legacy alias
                "beta_corr": r,
                "pearson_r": metrics["pearson_r"],
                "spearman_r": metrics["spearman_r"],
                "sign_concordance": metrics["sign_concordance"],
                "regression_slope": metrics["regression_slope"],
                "standardized_beta_rmse": metrics["standardized_beta_rmse"],
                "beta_corr_lo": lo,
                "beta_corr_hi": hi,
                "filter_frozen_on": "PanUKB_only",
                "filter": filter_name,
                "filter_training_columns": ",".join(train_cols),
                "n_after_filter": int(keep.sum()),
                "beta_corr_after_filter": r2,
                "beta_corr_after_lo": lo2,
                "beta_corr_after_hi": hi2,
                "delta_corr": (r2 - r) if np.isfinite(r2) and np.isfinite(r) else np.nan,
                "beta_corr_matched_random": rr,
                "delta_corr_random": (rr - r) if np.isfinite(rr) and np.isfinite(r) else np.nan,
                "status": status,
                "lift_status": "hg19ToHg38",
                **{f"qc_{k}": v for k, v in qc.items()},
                "note": "Filter frozen on Pan-UKB; allele QC ladder applied after liftover",
            }
        ]
    )
    out.to_csv(args.out, index=False, float_format="%.6g")
    pd.DataFrame([qc]).to_csv(TABLES / "external_page_qc_counts.csv", index=False)
    internal = TABLES / "internal_panukbb_concordance_sensitivity.csv"
    if internal.exists():
        comb = pd.concat([pd.read_csv(internal), out], ignore_index=True)
        comb.to_csv(TABLES / "external_sumstat_validation.csv", index=False, float_format="%.6g")
    print(out.to_string(index=False))
    print(f"Saved {args.out} status={status} n_variants={n}")


if __name__ == "__main__":
    main()
