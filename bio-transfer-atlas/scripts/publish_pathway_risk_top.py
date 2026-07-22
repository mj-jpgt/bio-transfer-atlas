"""
Build genome-wide pathway_risk_top summary from concatenated pathway parquet.
Marks FDR-significant pathways when enrichment stats are available.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TRAITS = ["CAD", "T2D", "BMI", "LDL"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Publish pathway_risk_top from genome-wide table.")
    p.add_argument(
        "--pathway-table",
        default=str(ROOT / "results/tables/pathway_risk_table_genomewide_genomewide.parquet"),
    )
    p.add_argument(
        "--enrichment-fdr",
        default=str(ROOT / "results/tables/pathway_enrichment_fdr_genomewide.csv"),
        help="Optional FDR table from pathway_enrichment_stats.py",
    )
    p.add_argument(
        "--out-csv",
        default=str(ROOT / "results/tables/pathway_risk_top_genomewide.csv"),
    )
    p.add_argument(
        "--out-summary",
        default=str(ROOT / "results/tables/pathway_risk_summary_genomewide.txt"),
    )
    p.add_argument("--top-n", type=int, default=8)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    path = Path(args.pathway_table)
    if not path.exists():
        raise SystemExit(f"Missing {path}")

    pr = pd.read_parquet(path)
    fdr_path = Path(args.enrichment_fdr)
    fdr = None
    if fdr_path.exists():
        fdr = pd.read_csv(fdr_path)
        keep = ["trait", "rid", "bh_q", "ld_block_perm_q", "fdr_significant", "fisher_p"]
        keep = [c for c in keep if c in fdr.columns]
        fdr = fdr[keep].drop_duplicates(["trait", "rid"])
        pr = pr.merge(fdr, on=["trait", "rid"], how="left")
        if "fdr_significant" in pr.columns:
            pr["fdr_significant"] = pr["fdr_significant"].fillna(False).astype(bool)
        else:
            pr["fdr_significant"] = False
        pr["claim_tier"] = pr["fdr_significant"].map(
            {True: "fdr_significant", False: "hypothesis_generating"}
        )
    else:
        pr["fdr_significant"] = False
        pr["claim_tier"] = "hypothesis_generating"

    lines = ["=" * 70, "TOP CROSS-ANCESTRY-UNSTABLE PATHWAYS (genome-wide)", "=" * 70]
    if fdr is not None:
        lines.append("(FDR table present: claim_tier = fdr_significant | hypothesis_generating)")
    else:
        lines.append("(No FDR table yet — all tops marked hypothesis_generating)")
    top_all = []

    for trait in TRAITS:
        sub = pr[(pr["trait"] == trait) & (pr["n_assoc"] >= 3)].copy()
        sub = sub.sort_values("mean_I2_assoc", ascending=False).head(args.top_n)
        lines.append(f"\n[{trait}]")
        if sub.empty:
            lines.append("  (no pathway with >=3 associated variants)")
            continue
        top_all.append(sub)
        for x in sub.itertuples(index=False):
            pathway = str(x.pathway)[:60]
            tier = getattr(x, "claim_tier", "hypothesis_generating")
            q = getattr(x, "bh_q", float("nan"))
            q_s = f" q={q:.3g}" if pd.notna(q) else ""
            lines.append(
                f"  [{tier}] I2={x.mean_I2_assoc:.3f}  n_assoc={int(x.n_assoc):3d}  "
                f"PBS={x.mean_PBS:.3f}{q_s}  {pathway}"
            )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    if top_all:
        pd.concat(top_all, ignore_index=True).to_csv(out_csv, index=False, float_format="%.4f")
    else:
        pd.DataFrame().to_csv(out_csv, index=False)

    Path(args.out_summary).write_text("\n".join(lines), encoding="utf-8")
    print(f"Saved {out_csv} ({sum(len(t) for t in top_all)} rows)")
    print(f"Saved {args.out_summary}")


if __name__ == "__main__":
    main()
