"""
FAIRGEN-Open Stage 6: Population Score Shifts
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

root = Path(__file__).resolve().parents[1]
DEFAULT_MATRIX = root / "data/processed/scores_grch38/score_matrix_grch38_genomewide.parquet"
out_dir = root / "results/tables"
out_dir.mkdir(parents=True, exist_ok=True)

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
PGS_TRAITS = {
    "PGS000018": "T2D",
    "PGS004696": "CAD",
    "PGS004698": "CAD",
    "PGS003897": "BMI",
    "PGS002853": "LDL",
    "PGS002858": "LDL",
    "PGS003092": "BMI",
    "PGS000014": "CAD",
    "PGS004840": "T2D",
    "PGS000191": "WBC",
    "PGS004133": "RA",
    "PGS001288": "IBD",
}
PGS_IDS = list(PGS_TRAITS.keys())


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute population score shifts from score matrix.")
    p.add_argument("--score-matrix", default=str(DEFAULT_MATRIX))
    p.add_argument("--tag", default="genomewide", help="Suffix for output tables")
    return p.parse_args()


def kendall_w(ranks: np.ndarray) -> float:
    n_raters, n_subjects = ranks.shape
    r = ranks.sum(axis=0)
    r_bar = r.mean()
    s = np.sum((r - r_bar) ** 2)
    return float(12 * s / (n_raters**2 * (n_subjects**3 - n_subjects)))


def main() -> None:
    args = parse_args()
    tag = args.tag
    score_path = Path(args.score_matrix)
    if not score_path.exists():
        score_path = root / "data/processed/scores_grch38/score_matrix_grch38.parquet"
    print(f"Loading score matrix: {score_path}")
    scores = pd.read_parquet(score_path)
    scores = scores.dropna(subset=["super_pop"])
    print(f"  {len(scores)} samples with pop labels")

    sp_rows = []
    for pgs in PGS_IDS:
        if pgs not in scores.columns:
            continue
        for sp in SUPERPOPS:
            sub = scores[scores["super_pop"] == sp][pgs].dropna()
            if len(sub) < 10:
                continue
            sp_rows.append({
                "pgs_id": pgs,
                "trait": PGS_TRAITS[pgs],
                "super_pop": sp,
                "n": len(sub),
                "mean": sub.mean(),
                "sd": sub.std(),
                "median": sub.median(),
                "p10": sub.quantile(0.10),
                "p90": sub.quantile(0.90),
                "iqr": sub.quantile(0.75) - sub.quantile(0.25),
            })

    sp_df = pd.DataFrame(sp_rows)
    eur_stats = sp_df[sp_df["super_pop"] == "EUR"][["pgs_id", "mean", "sd"]].rename(
        columns={"mean": "mean_EUR", "sd": "sd_EUR"}
    )
    sp_df = sp_df.merge(eur_stats, on="pgs_id", how="left")
    sp_df["delta_EUR"] = (sp_df["mean"] - sp_df["mean_EUR"]) / sp_df["sd_EUR"]
    sp_df["rank_by_mean"] = sp_df.groupby("pgs_id")["mean"].rank(ascending=False).astype(int)
    sp_df.to_parquet(out_dir / f"score_shifts_superpop_{tag}.parquet", index=False)

    pop_rows = []
    for pgs in PGS_IDS:
        if pgs not in scores.columns:
            continue
        for pop in scores["pop"].dropna().unique():
            sub = scores[scores["pop"] == pop][pgs].dropna()
            if len(sub) < 5:
                continue
            pop_rows.append({
                "pgs_id": pgs,
                "pop": pop,
                "super_pop": scores[scores["pop"] == pop]["super_pop"].iloc[0],
                "n": len(sub),
                "mean": sub.mean(),
                "sd": sub.std(),
                "p10": sub.quantile(0.10),
                "p90": sub.quantile(0.90),
            })
    pop_df = pd.DataFrame(pop_rows)
    eur_pop = pop_df[pop_df["super_pop"] == "EUR"].groupby("pgs_id")["mean"].mean().rename("mean_EUR_mean")
    pop_df = pop_df.merge(eur_pop, on="pgs_id", how="left")
    pop_df["delta_EUR"] = (pop_df["mean"] - pop_df["mean_EUR_mean"]) / pop_df.merge(
        sp_df[sp_df["super_pop"] == "EUR"][["pgs_id", "sd_EUR"]], on="pgs_id", how="left"
    )["sd_EUR"]
    pop_df.to_parquet(out_dir / f"score_shifts_pop_{tag}.parquet", index=False)

    rank_matrix = sp_df.pivot(index="super_pop", columns="pgs_id", values="rank_by_mean")
    rank_matrix = rank_matrix.reindex(SUPERPOPS)
    rank_matrix.to_parquet(out_dir / f"rank_instability_{tag}.parquet")
    w = kendall_w(rank_matrix.values.astype(float))
    print(f"Kendall's W across PGS: {w:.4f}")

    summary = sp_df[
        ["pgs_id", "trait", "super_pop", "n", "mean", "sd", "delta_EUR", "rank_by_mean"]
    ].sort_values(["pgs_id", "rank_by_mean"])
    summary.to_csv(out_dir / f"score_shift_summary_{tag}.csv", index=False, float_format="%.5f")
    print(f"Saved score_shift_summary_{tag}.csv ({len(summary)} rows)")


if __name__ == "__main__":
    main()
