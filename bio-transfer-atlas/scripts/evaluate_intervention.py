"""
Phase 18.5: Evaluate genome-wide intervention impact vs baseline and negative controls.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))

from intervention_common import (
    GW_INT_ROOT,
    GW_MASTER,
    GW_PREDS,
    GW_SCORE_MATRIX,
    GW_TAG,
    INTERVENTION_MODES,
    PGS_IDS,
    PGS_TRAITS,
    PLINK2,
    ROOT,
    SUPERPOPS,
    available_score_chroms,
    load_genomewide_weights,
    pfile_for_chrom,
)

OUT_DIR = ROOT / "results/tables"
INT_SCORES = ROOT / "data/processed/scores_grch38_intervention_genomewide"
DEFAULT_CSV = OUT_DIR / f"intervention_results.{GW_TAG}.csv"
DEFAULT_SUMMARY = OUT_DIR / f"intervention_summary.{GW_TAG}.txt"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Evaluate genome-wide intervention outcomes.")
    p.add_argument("--baseline", default=str(GW_SCORE_MATRIX))
    p.add_argument("--intervention-scores", default=str(INT_SCORES))
    p.add_argument("--labels", default=str(GW_MASTER))
    p.add_argument("--intervention-root", default=str(GW_INT_ROOT))
    p.add_argument("--modes", default=",".join(INTERVENTION_MODES))
    p.add_argument("--tag", default=GW_TAG, help="Suffix matching score_matrix_{mode}_{tag}.parquet")
    p.add_argument("--n-boot", type=int, default=100, help="Bootstrap resamples for MAD CI")
    p.add_argument("--out-csv", default=str(DEFAULT_CSV))
    p.add_argument("--out-summary", default=str(DEFAULT_SUMMARY))
    return p.parse_args()


def superpop_stats(scores: pd.DataFrame, pgs_id: str) -> pd.DataFrame:
    rows = []
    for sp in SUPERPOPS:
        sub = scores.loc[scores["super_pop"] == sp, pgs_id].dropna()
        if len(sub) < 10:
            continue
        rows.append({"super_pop": sp, "mean": sub.mean(), "sd": sub.std()})
    sp_df = pd.DataFrame(rows)
    eur = sp_df.loc[sp_df["super_pop"] == "EUR"]
    if eur.empty:
        sp_df["delta_EUR"] = np.nan
        return sp_df
    mean_eur = float(eur["mean"].iloc[0])
    sd_eur = float(eur["sd"].iloc[0])
    sp_df["delta_EUR"] = (sp_df["mean"] - mean_eur) / sd_eur if sd_eur > 0 else np.nan
    return sp_df


def mean_abs_delta_eur(sp_df: pd.DataFrame) -> float:
    non_eur = sp_df[sp_df["super_pop"] != "EUR"]
    if non_eur.empty:
        return float("nan")
    return float(non_eur["delta_EUR"].abs().mean())


def bootstrap_mad(scores: pd.DataFrame, pgs_id: str, n_boot: int = 200, seed: int = 719) -> tuple[float, float, float]:
    """Bootstrap CI for mean |delta_EUR| over individuals (resample within superpops)."""
    rng = np.random.default_rng(seed)
    base = mean_abs_delta_eur(superpop_stats(scores, pgs_id))
    vals = []
    by_sp = {sp: scores.loc[scores["super_pop"] == sp, pgs_id].dropna().to_numpy() for sp in SUPERPOPS}
    for _ in range(n_boot):
        rows = []
        for sp, arr in by_sp.items():
            if len(arr) < 10:
                continue
            idx = rng.integers(0, len(arr), len(arr))
            rows.append(pd.DataFrame({"super_pop": sp, pgs_id: arr[idx]}))
        if not rows:
            continue
        boot = pd.concat(rows, ignore_index=True)
        vals.append(mean_abs_delta_eur(superpop_stats(boot, pgs_id)))
    if not vals:
        return base, float("nan"), float("nan")
    return base, float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def score_tail_ancestry(scores: pd.DataFrame, pgs_id: str, frac: float = 0.05) -> dict[str, float]:
    """Fraction of each superpop among top |score| percentile (absolute extreme)."""
    s = scores[["super_pop", pgs_id]].dropna()
    if s.empty:
        return {sp: float("nan") for sp in SUPERPOPS}
    thr = s[pgs_id].abs().quantile(1.0 - frac)
    tail = s[s[pgs_id].abs() >= thr]
    if tail.empty:
        return {sp: float("nan") for sp in SUPERPOPS}
    counts = tail["super_pop"].value_counts(normalize=True)
    return {sp: float(counts.get(sp, 0.0)) for sp in SUPERPOPS}


def kendall_w(rank_matrix: np.ndarray) -> float:
    n_raters, n_subjects = rank_matrix.shape
    if n_subjects < 2:
        return float("nan")
    r = rank_matrix.sum(axis=0)
    r_bar = r.mean()
    s = np.sum((r - r_bar) ** 2)
    return float(12 * s / (n_raters**2 * (n_subjects**3 - n_subjects)))


def rank_instability(scores: pd.DataFrame, pgs_ids: list[str]) -> float:
    ranks = []
    for pgs in pgs_ids:
        if pgs not in scores.columns:
            continue
        sp_df = superpop_stats(scores, pgs)
        sp_df = sp_df.set_index("super_pop").reindex(SUPERPOPS)
        sp_df["rank"] = sp_df["mean"].rank(ascending=False)
        ranks.append(sp_df["rank"].to_numpy(dtype=float))
    if len(ranks) < 2:
        return float("nan")
    return kendall_w(np.vstack(ranks))


def retained_variant_ids(pgs_id: str, mode: str, int_root: Path, cache: dict) -> set[str]:
    key = (pgs_id, mode)
    if key in cache:
        return cache[key]
    path = int_root / pgs_id / f"{mode}.tsv"
    if not path.exists():
        cache[key] = set()
        return cache[key]
    df = pd.read_csv(path, sep="\t", usecols=["variant_id"], dtype=str)
    cache[key] = set(df["variant_id"])
    return cache[key]


def load_labels(labels_path: Path) -> pd.DataFrame:
    """Stream only needed columns; prefer associated rows when available."""
    import pyarrow.dataset as ds

    cols = ["variant_id", "trait", "I2", "sign_concordance"]
    dataset = ds.dataset(str(labels_path), format="parquet")
    available = set(dataset.schema.names)
    use_cols = [c for c in cols if c in available]
    filt = None
    if "associated" in available:
        filt = ds.field("associated") == True  # noqa: E712
        print("  loading associated label rows only ...", flush=True)
    scanner = dataset.scanner(columns=use_cols, filter=filt, batch_size=1_000_000)
    chunks = []
    n = 0
    for batch in scanner.to_batches():
        chunks.append(batch.to_pandas().drop_duplicates(["variant_id", "trait"]))
        n += len(chunks[-1])
        if n % 2_000_000 < 1_000_000:
            print(f"  labels loaded {n:,} ...", flush=True)
    if not chunks:
        return pd.DataFrame(columns=cols)
    out = pd.concat(chunks, ignore_index=True).drop_duplicates(["variant_id", "trait"])
    print(f"  labels: {len(out):,} variant×trait rows", flush=True)
    return out


def concordance_on_retained(
    labels_by_trait: dict[str, pd.DataFrame],
    trait: str,
    retained: set[str],
) -> tuple[float, float, int]:
    sub = labels_by_trait.get(trait)
    if sub is None or sub.empty or not retained:
        return float("nan"), float("nan"), 0
    # Boolean mask via map is much faster than isin on multi-million sets
    mask = sub["variant_id"].map(retained.__contains__)
    hit = sub.loc[mask]
    if hit.empty:
        return float("nan"), float("nan"), 0
    return float(hit["I2"].median()), float(hit["sign_concordance"].mean()), len(hit)


def main() -> None:
    args = parse_args()
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    out_dir = Path(args.out_csv).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading baseline scores ...", flush=True)
    baseline = pd.read_parquet(args.baseline).dropna(subset=["super_pop"])
    print(f"Loading labels from {args.labels} ...", flush=True)
    labels = load_labels(Path(args.labels))
    labels_by_trait = {t: g.reset_index(drop=True) for t, g in labels.groupby("trait", sort=False)}
    int_root = Path(args.intervention_root)
    int_scores_root = Path(args.intervention_scores)
    retained_cache: dict = {}
    weights_cache: dict[str, pd.DataFrame] = {}

    def weights_for(pgs: str) -> pd.DataFrame:
        if pgs not in weights_cache:
            print(f"  caching weights {pgs} ...", flush=True)
            weights_cache[pgs] = load_genomewide_weights(pgs)
        return weights_cache[pgs]

    baseline_metrics = {}
    baseline_conc = {}
    baseline_tail = {}
    print("Computing baseline metrics ...", flush=True)
    for pgs in PGS_IDS:
        if pgs not in baseline.columns:
            continue
        print(f"  baseline {pgs}", flush=True)
        mad, mad_lo, mad_hi = bootstrap_mad(baseline, pgs, n_boot=args.n_boot)
        baseline_metrics[pgs] = {
            "mean_abs_delta_EUR": mad,
            "mean_abs_delta_EUR_lo": mad_lo,
            "mean_abs_delta_EUR_hi": mad_hi,
            "kendall_w": rank_instability(baseline, PGS_IDS),
        }
        orig_ids = set(weights_for(pgs)["variant_id"].astype(str))
        med_i2_b, mean_sc_b, n_lab_b = concordance_on_retained(
            labels_by_trait, PGS_TRAITS[pgs], orig_ids
        )
        baseline_conc[pgs] = {
            "median_I2": med_i2_b,
            "mean_sign_concordance": mean_sc_b,
            "n_labeled": n_lab_b,
        }
        baseline_tail[pgs] = score_tail_ancestry(baseline, pgs)

    rows = []
    for mode in modes:
        score_path = int_scores_root / f"score_matrix_{mode}_{args.tag}.parquet"
        if not score_path.exists():
            print(f"skip mode {mode}: missing {score_path}", flush=True)
            continue
        print(f"Evaluating mode {mode} ...", flush=True)
        scores = pd.read_parquet(score_path).dropna(subset=["super_pop"])
        w_instab = rank_instability(scores, PGS_IDS)

        for pgs in PGS_IDS:
            if pgs not in scores.columns or pgs not in baseline_metrics:
                continue
            trait = PGS_TRAITS[pgs]
            mad, mad_lo, mad_hi = bootstrap_mad(scores, pgs, n_boot=args.n_boot)
            base_mad = baseline_metrics[pgs]["mean_abs_delta_EUR"]
            retained = retained_variant_ids(pgs, mode, int_root, retained_cache)
            med_i2, mean_sc, n_lab = concordance_on_retained(labels_by_trait, trait, retained)
            base_i2 = baseline_conc[pgs]["median_I2"]
            base_sc = baseline_conc[pgs]["mean_sign_concordance"]

            orig_w = weights_for(pgs)
            wsum_orig = orig_w["effect_weight"].abs().sum()
            if mode.startswith("reweight"):
                mod = pd.read_csv(int_root / pgs / f"{mode}.tsv", sep="\t")
                w_ret = float(mod["effect_weight"].abs().sum() / wsum_orig) if wsum_orig else np.nan
            elif mode == "flag":
                w_ret = 1.0
            else:
                kept = orig_w[orig_w["variant_id"].isin(retained)]
                w_ret = float(kept["effect_weight"].abs().sum() / wsum_orig) if wsum_orig else np.nan

            rows.append({
                "pgs_id": pgs, "trait": trait, "mode": mode,
                "metric": "mean_abs_delta_EUR", "value": mad,
                "value_lo": mad_lo, "value_hi": mad_hi,
                "baseline_value": base_mad,
                "baseline_value_lo": baseline_metrics[pgs]["mean_abs_delta_EUR_lo"],
                "baseline_value_hi": baseline_metrics[pgs]["mean_abs_delta_EUR_hi"],
                "delta_vs_baseline": mad - base_mad if pd.notna(mad) and pd.notna(base_mad) else np.nan,
                "reduction": base_mad - mad if pd.notna(mad) and pd.notna(base_mad) else np.nan,
            })
            rows.append({
                "pgs_id": pgs, "trait": trait, "mode": mode,
                "metric": "kendall_w", "value": w_instab,
                "baseline_value": baseline_metrics[pgs]["kendall_w"],
                "delta_vs_baseline": w_instab - baseline_metrics[pgs]["kendall_w"],
                "reduction": baseline_metrics[pgs]["kendall_w"] - w_instab,
            })
            rows.append({
                "pgs_id": pgs, "trait": trait, "mode": mode,
                "metric": "median_I2_retained", "value": med_i2,
                "baseline_value": base_i2,
                "delta_vs_baseline": med_i2 - base_i2 if pd.notna(med_i2) and pd.notna(base_i2) else np.nan,
                "reduction": base_i2 - med_i2 if pd.notna(med_i2) and pd.notna(base_i2) else np.nan,
                "n_labeled_retained": n_lab,
            })
            rows.append({
                "pgs_id": pgs, "trait": trait, "mode": mode,
                "metric": "mean_sign_concordance_retained", "value": mean_sc,
                "baseline_value": base_sc,
                "delta_vs_baseline": mean_sc - base_sc if pd.notna(mean_sc) and pd.notna(base_sc) else np.nan,
                "reduction": mean_sc - base_sc if pd.notna(mean_sc) and pd.notna(base_sc) else np.nan,
                "n_labeled_retained": n_lab,
            })
            rows.append({
                "pgs_id": pgs, "trait": trait, "mode": mode,
                "metric": "weight_retained_frac", "value": w_ret,
                "baseline_value": 1.0,
                "delta_vs_baseline": w_ret - 1.0 if pd.notna(w_ret) else np.nan,
                "reduction": 1.0 - w_ret if pd.notna(w_ret) else np.nan,
            })
            n_ret = (
                float(len(retained) / max(len(orig_w), 1))
                if isinstance(retained, set)
                else np.nan
            )
            rows.append({
                "pgs_id": pgs, "trait": trait, "mode": mode,
                "metric": "variant_frac_retained", "value": n_ret,
                "baseline_value": 1.0,
                "claim": "ancestry_mean_separation_not_portability",
            })
            # Within-ancestry variance retained + corr(edited, original)
            for sp in SUPERPOPS:
                b = baseline.loc[baseline["super_pop"] == sp, ["sample_id", pgs]].rename(
                    columns={pgs: "base"}
                )
                e = scores.loc[scores["super_pop"] == sp, ["sample_id", pgs]].rename(
                    columns={pgs: "edit"}
                )
                m = b.merge(e, on="sample_id")
                if len(m) < 10:
                    continue
                vb, ve = float(m["base"].var()), float(m["edit"].var())
                rows.append({
                    "pgs_id": pgs, "trait": trait, "mode": mode,
                    "metric": f"score_var_retained_{sp}",
                    "value": (ve / vb) if vb > 0 else np.nan,
                    "claim": "ancestry_mean_separation_not_portability",
                })
                rows.append({
                    "pgs_id": pgs, "trait": trait, "mode": mode,
                    "metric": f"corr_edited_original_{sp}",
                    "value": float(np.corrcoef(m["base"], m["edit"])[0, 1]),
                    "claim": "ancestry_mean_separation_not_portability",
                })

            # Score-tail ancestry composition (top 5% |score|)
            tail = score_tail_ancestry(scores, pgs)
            base_tail = baseline_tail[pgs]
            for sp in SUPERPOPS:
                rows.append({
                    "pgs_id": pgs, "trait": trait, "mode": mode,
                    "metric": f"tail5_frac_{sp}",
                    "value": tail[sp],
                    "baseline_value": base_tail[sp],
                    "delta_vs_baseline": (
                        tail[sp] - base_tail[sp]
                        if pd.notna(tail[sp]) and pd.notna(base_tail[sp])
                        else np.nan
                    ),
                    "reduction": np.nan,
                })
            # L1 distance of tail composition vs baseline
            l1 = sum(
                abs(tail[sp] - base_tail[sp])
                for sp in SUPERPOPS
                if pd.notna(tail[sp]) and pd.notna(base_tail[sp])
            )
            rows.append({
                "pgs_id": pgs, "trait": trait, "mode": mode,
                "metric": "tail5_ancestry_L1",
                "value": l1,
                "baseline_value": 0.0,
                "delta_vs_baseline": l1,
                "reduction": np.nan,
            })
            print(f"  {mode}/{pgs} MAD={mad:.4f} reduction={base_mad - mad:.4f}", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(args.out_csv, index=False)
    print(f"Saved -> {args.out_csv}", flush=True)

    # Matched random Monte Carlo + LOSO (AF-expected MAD; scored MAD when matrices exist)
    try:
        from run_intervention_matched_random_controls import main as matched_main
        import sys as _sys

        print("Running matched random controls + LOSO ...", flush=True)
        _argv = _sys.argv
        _sys.argv = ["run_intervention_matched_random_controls.py"]
        try:
            matched_main()
        finally:
            _sys.argv = _argv
    except Exception as exc:
        print(f"matched controls skip: {exc}", flush=True)

    lines = [f"Phase 18 {args.tag} intervention summary", "=" * 40]
    mad = res[res["metric"] == "mean_abs_delta_EUR"].copy()
    for mode in ["reweight_linear", "filter_10", "random"]:
        sub = mad[mad["mode"] == mode]
        if sub.empty:
            continue
        lines.append(f"\n{mode}:")
        for _, r in sub.iterrows():
            lines.append(
                f"  {r['pgs_id']}: MAD={r['value']:.4f} "
                f"(baseline {r['baseline_value']:.4f}, reduction {r['reduction']:.4f})"
            )

    summary_text = "\n".join(lines)
    Path(args.out_summary).write_text(summary_text)
    print(summary_text, flush=True)
    print(f"\nSaved -> {args.out_summary}", flush=True)


if __name__ == "__main__":
    main()
