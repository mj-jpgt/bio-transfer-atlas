"""
FAIRGEN-Open Stage 7 (cont.): Genetic-distance sensitivity
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

root = Path(__file__).resolve().parents[1]
proc = root / "data/processed"
out_dir = root / "results/tables"
DEFAULT_MATRIX = proc / "scores_grch38/score_matrix_grch38_genomewide.parquet"

N_PC = 10
N_BOOT = 1000
SEED = 719

PGS_TRAITS = {
    "PGS000018": "T2D", "PGS004840": "T2D",
    "PGS004696": "CAD", "PGS004698": "CAD", "PGS000014": "CAD",
    "PGS003897": "BMI", "PGS003092": "BMI",
    "PGS002853": "LDL", "PGS002858": "LDL",
    "PGS000191": "WBC", "PGS004133": "RA", "PGS001288": "IBD",
}
PGS_IDS = list(PGS_TRAITS)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Distance sensitivity from score matrix.")
    p.add_argument("--score-matrix", default=str(DEFAULT_MATRIX))
    p.add_argument("--tag", default="genomewide")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    tag = args.tag
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    score_path = Path(args.score_matrix)
    if not score_path.exists():
        score_path = proc / "scores_grch38/score_matrix_grch38.parquet"
    scores = pd.read_parquet(score_path)
    scores = scores.dropna(subset=["super_pop", "pop"])
    pca = pd.read_parquet(proc / "ancestry_pcs_genome_wide.parquet")
    pc_cols = [f"PC{i}" for i in range(1, N_PC + 1)]

    df = scores.merge(pca[["sample_id"] + pc_cols], on="sample_id", how="inner")
    centroids = df.groupby("pop")[pc_cols].mean()
    pop_super = df.groupby("pop")["super_pop"].first()
    eur_centroid = df[df["super_pop"] == "EUR"][pc_cols].mean().values

    dist = np.sqrt(((centroids[pc_cols].values - eur_centroid) ** 2).sum(axis=1))
    dist_df = pd.DataFrame({
        "pop": centroids.index,
        "super_pop": pop_super.reindex(centroids.index).values,
        "dist_to_EUR": dist,
    }).sort_values("dist_to_EUR")
    dist_df.to_parquet(out_dir / f"pop_genetic_distance_{tag}.parquet", index=False)

    def pop_shift_table(pgs: str) -> pd.DataFrame:
        eur_vals = df[df["super_pop"] == "EUR"][pgs].dropna()
        mu_eur, sd_eur = eur_vals.mean(), eur_vals.std()
        rows = []
        for pop in centroids.index:
            v = df[df["pop"] == pop][pgs].dropna()
            if len(v) < 5:
                continue
            rows.append({"pop": pop, "shift": (v.mean() - mu_eur) / sd_eur})
        return pd.DataFrame(rows).merge(dist_df, on="pop", how="left")

    def fit_slope(d: np.ndarray, y: np.ndarray) -> tuple[float, float]:
        res = stats.linregress(d, y)
        return res.slope, res.rvalue ** 2

    results = []
    for pgs in PGS_IDS:
        if pgs not in df.columns:
            continue
        t = pop_shift_table(pgs)
        t = t[t["super_pop"] != "EUR"]
        d = t["dist_to_EUR"].values
        y = t["shift"].abs().values
        if len(d) < 4:
            continue
        slope, r2 = fit_slope(d, y)
        boots = []
        n = len(d)
        for _ in range(N_BOOT):
            idx = rng.integers(0, n, n)
            boots.append(fit_slope(d[idx], y[idx])[0])
        lo, hi = np.percentile(boots, [2.5, 97.5])
        results.append({
            "pgs_id": pgs, "trait": PGS_TRAITS[pgs],
            "slope": slope, "r2": r2, "slope_lo": lo, "slope_hi": hi, "n_pops": n,
        })

    res_df = pd.DataFrame(results).sort_values("slope", ascending=False)
    res_df.to_parquet(out_dir / f"distance_sensitivity_{tag}.parquet", index=False)
    res_df.to_csv(out_dir / f"distance_sensitivity_summary_{tag}.csv", index=False, float_format="%.5f")
    print(f"Saved distance_sensitivity_{tag} ({len(res_df)} PGS)")


if __name__ == "__main__":
    main()
