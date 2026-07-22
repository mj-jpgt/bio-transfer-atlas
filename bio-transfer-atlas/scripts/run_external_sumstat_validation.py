"""
Phase C2: Open external validation scaffolding (GBMI / PAGE sumstat concordance).

Downloads are optional; if local sumstats exist, computes cross-ancestry sign concordance
and effect-size correlation vs Pan-UKB EUR, and intervention reweight impact on concordance.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--external",
        default=str(ROOT / "data/raw/external_sumstats"),
        help="Dir with {source}_{trait}_{anc}.parquet columns: variant_id,beta,se,pval",
    )
    p.add_argument(
        "--predictions",
        default=str(ROOT / "data/modeling/variant_portability_predictions.genomewide.parquet"),
    )
    p.add_argument(
        "--out",
        default=str(ROOT / "results/tables/external_sumstat_validation.csv"),
    )
    p.add_argument("--top-risk-frac", type=float, default=0.1)
    return p.parse_args()


def bootstrap_corr(x: np.ndarray, y: np.ndarray, n: int = 200, seed: int = 719) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 30:
        return float("nan"), float("nan"), float("nan")
    base = float(np.corrcoef(x, y)[0, 1])
    vals = []
    for _ in range(n):
        idx = rng.integers(0, len(x), len(x))
        vals.append(float(np.corrcoef(x[idx], y[idx])[0, 1]))
    return base, float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main() -> None:
    args = parse_args()
    ext_dir = Path(args.external)
    ext_dir.mkdir(parents=True, exist_ok=True)

    # Manifest of expected open resources
    manifest = [
        {
            "source": "GBMI",
            "note": "Place GBMI multi-ancestry sumstats as gbmi_{trait}_{anc}.parquet",
            "url": "https://www.globalbiobankmeta.org/",
        },
        {
            "source": "PAGE_Wojcik2019",
            "note": "GWAS Catalog PAGE sumstats as page_{trait}_{anc}.parquet",
            "url": "https://www.ebi.ac.uk/gwas/",
        },
        {
            "source": "AllOfUs",
            "note": "Parallel track via Researcher Workbench — not required for open path",
            "url": "https://workbench.researchallofus.org/",
        },
    ]
    (ext_dir / "README_external_sumstats.md").write_text(
        "# External sumstats for validation\n\n"
        + "\n".join(f"- **{m['source']}**: {m['note']} ({m['url']})" for m in manifest)
        + "\n\nRequired columns: variant_id, beta, se, pval (GRCh38 preferred).\n",
        encoding="utf-8",
    )

    files = list(ext_dir.glob("*.parquet"))
    rows = []
    if not files:
        rows.append(
            {
                "status": "awaiting_sumstats",
                "n_files": 0,
                "message": "No external parquet files yet; README written. Open path ready.",
            }
        )
        out = pd.DataFrame(rows)
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out, index=False)
        print(out.to_string(index=False))
        print(f"Saved {args.out}")
        return

    # Load risk scores
    import pyarrow.dataset as ds

    dataset = ds.dataset(str(args.predictions), format="parquet")
    names = set(dataset.schema.names)
    risk_col = next(
        (
            c
            for c in ["y_prob", "pred_prob", "prob_high_I2", "portability_risk"]
            if c in names
        ),
        None,
    )
    if risk_col is None:
        raise SystemExit("No risk column in predictions")
    chunks = []
    for batch in dataset.scanner(columns=["variant_id", risk_col], batch_size=500_000).to_batches():
        chunks.append(batch.to_pandas())
    risk = pd.concat(chunks, ignore_index=True).groupby("variant_id", as_index=False)[risk_col].max()
    thr = risk[risk_col].quantile(1 - args.top_risk_frac)
    high_risk = set(risk.loc[risk[risk_col] >= thr, "variant_id"])

    # Pair EUR vs non-EUR files by trait
    by_key: dict[str, dict[str, Path]] = {}
    for f in files:
        # {source}_{trait}_{anc}.parquet
        parts = f.stem.split("_")
        if len(parts) < 3:
            continue
        source, trait, anc = parts[0], parts[1], parts[2]
        by_key.setdefault(f"{source}_{trait}", {})[anc.upper()] = f

    for key, ancs in by_key.items():
        if "EUR" not in ancs:
            continue
        eur = pd.read_parquet(ancs["EUR"])
        for anc, path in ancs.items():
            if anc == "EUR":
                continue
            other = pd.read_parquet(path)
            m = eur.merge(other, on="variant_id", suffixes=("_eur", "_anc"))
            if m.empty:
                continue
            r, lo, hi = bootstrap_corr(m["beta_eur"].to_numpy(), m["beta_anc"].to_numpy())
            # After dropping top-risk variants
            m2 = m[~m["variant_id"].isin(high_risk)]
            r2, lo2, hi2 = bootstrap_corr(m2["beta_eur"].to_numpy(), m2["beta_anc"].to_numpy())
            sign_conc = float(np.mean(np.sign(m["beta_eur"]) == np.sign(m["beta_anc"])))
            sign_conc2 = float(np.mean(np.sign(m2["beta_eur"]) == np.sign(m2["beta_anc"])))
            rows.append(
                {
                    "pair": key,
                    "anc": anc,
                    "n": len(m),
                    "beta_corr": r,
                    "beta_corr_lo": lo,
                    "beta_corr_hi": hi,
                    "beta_corr_after_filter": r2,
                    "beta_corr_after_lo": lo2,
                    "beta_corr_after_hi": hi2,
                    "delta_corr": r2 - r if pd.notna(r2) and pd.notna(r) else np.nan,
                    "sign_concordance": sign_conc,
                    "sign_concordance_after_filter": sign_conc2,
                    "status": "ok",
                }
            )

    out = pd.DataFrame(rows) if rows else pd.DataFrame([{"status": "no_eur_pairs"}])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, float_format="%.6g")
    print(out.to_string(index=False))
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
