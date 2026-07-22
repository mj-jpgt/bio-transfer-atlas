#!/usr/bin/env python3
"""
Duffy positive control with allele audit, homozygous contrast, and ACKR1 score split.

Primary AFR biology contrast: duffy_null_dose==2 vs dose<2 (with bootstrap CIs).
Dominant null_gt0 / null_gt_ge1 kept as complementary labels (not the same group).
ACKR1 decomposition: score_snp + score_window_rest + score_outside ≈ score_total.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from bta_plink import plink2_bin  # noqa: E402
from intervention_common import (  # noqa: E402
    DUFFY_CHR,
    DUFFY_NULL_ALLELE,
    DUFFY_POS,
    DUFFY_RSID,
    DUFFY_WINDOW_BP,
    load_genomewide_weights,
)

WBC_PGS = "PGS000191"
PANEL = ROOT / "data/processed/sample_metadata_grch38.parquet"
PFILE_ROOT = ROOT / "data/interim/1000g_grch38"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--score-matrix",
        default=str(ROOT / "data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet"),
    )
    p.add_argument(
        "--predictions",
        default=str(ROOT / "data/modeling/variant_portability_predictions.genomewide.parquet"),
    )
    p.add_argument("--pfile", default=str(ROOT / "data/interim/1000g_grch38/chr1.score"))
    p.add_argument("--out-dir", default=str(ROOT / "results/tables"))
    p.add_argument("--fig-dir", default=str(ROOT / "results/figures"))
    return p.parse_args()


def extract_duffy_dosages(pfile: Path, work: Path) -> tuple[pd.DataFrame, dict]:
    work.mkdir(parents=True, exist_ok=True)
    out = work / "duffy_window"
    cmd = [
        plink2_bin(),
        "--pfile",
        str(pfile),
        "--chr",
        DUFFY_CHR,
        "--from-bp",
        str(DUFFY_POS - 50_000),
        "--to-bp",
        str(DUFFY_POS + 50_000),
        "--export",
        "A",
        "--out",
        str(out),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    raw = Path(str(out) + ".raw")
    if r.returncode != 0 or not raw.exists():
        raise SystemExit(f"plink export failed: {r.stderr[-500:]}")
    df = pd.read_csv(raw, sep=r"\s+")
    id_col = "IID" if "IID" in df.columns else "#IID"
    dose_cols = [c for c in df.columns if c not in ("FID", "IID", "#IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
    prefer = [c for c in dose_cols if str(DUFFY_POS) in c]
    col = prefer[0] if prefer else dose_cols[0]
    counted = col.rsplit("_", 1)[-1].upper() if "_" in col else "?"
    id_part = col.rsplit("_", 1)[0]
    parts = id_part.split(":")
    ref = parts[2].upper() if len(parts) >= 4 else "?"
    alt = parts[3].upper() if len(parts) >= 4 else "?"
    dose = pd.to_numeric(df[col], errors="coerce")
    if counted == DUFFY_NULL_ALLELE.upper():
        null_dose = dose
        flipped = False
    elif counted in (ref, alt) and DUFFY_NULL_ALLELE.upper() in (ref, alt):
        null_dose = 2.0 - dose
        flipped = True
    else:
        null_dose = dose
        flipped = False
    out_df = pd.DataFrame(
        {
            "sample_id": df[id_col].astype(str),
            "plink_counted_allele_dose": dose,
            "duffy_null_dose": null_dose,
            "duffy_snp_col": col,
        }
    )
    out_df["duffy_null_homozygous"] = (out_df["duffy_null_dose"].round() == 2).astype(int)
    audit = {
        "rsid": DUFFY_RSID,
        "chrom": DUFFY_CHR,
        "pos_grch38": DUFFY_POS,
        "ref": ref,
        "alt": alt,
        "plink_column": col,
        "plink_counted_allele": counted,
        "duffy_null_allele_target": DUFFY_NULL_ALLELE,
        "dosage_flipped_to_null": flipped,
        "interpretation": (
            "duffy_null_dose = copies of Duffy-null allele "
            f"({DUFFY_NULL_ALLELE}); expect high in AFR, low in EUR"
        ),
    }
    return out_df, audit


def _boot_mean_ci(vals: np.ndarray, n_boot: int = 400, seed: int = 719) -> tuple[float, float, float]:
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return float("nan"), float("nan"), float("nan")
    mean = float(np.mean(vals))
    if len(vals) == 1:
        return mean, mean, mean
    rng = np.random.default_rng(seed)
    boots = [float(np.mean(vals[rng.integers(0, len(vals), len(vals))])) for _ in range(n_boot)]
    return mean, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def genotype_strata(scores: pd.DataFrame, dose: pd.DataFrame) -> pd.DataFrame:
    m = scores.merge(dose, on="sample_id", how="inner")
    m = m.dropna(subset=["duffy_null_dose", "super_pop", WBC_PGS])
    m["duffy_null_gt"] = m["duffy_null_dose"].round().clip(0, 2).astype(int)
    rows = []
    for sp in ["AFR", "EUR", "AMR", "EAS", "SAS"]:
        sub = m[m["super_pop"] == sp]
        if sub.empty:
            continue
        for gt in [0, 1, 2]:
            g = sub[sub["duffy_null_gt"] == gt]
            n = len(g)
            if n == 0:
                continue
            vals = g[WBC_PGS].to_numpy(float)
            mean, lo, hi = _boot_mean_ci(vals)
            rows.append(
                {
                    "super_pop": sp,
                    "duffy_null_gt": gt,
                    "contrast": "genotype_dose",
                    "n": n,
                    "mean_wbc_pgs": mean,
                    "sd_wbc_pgs": float(np.std(vals, ddof=1)) if n > 1 else 0.0,
                    "mean_lo": lo,
                    "mean_hi": hi,
                    "underpowered": n < 10,
                    "allele_definition": f"copies of {DUFFY_NULL_ALLELE} (Duffy-null)",
                }
            )
        # Dominant complements (not the same group)
        for label, mask in [
            ("null_gt0", sub["duffy_null_gt"] == 0),
            ("null_gt_ge1", sub["duffy_null_gt"] >= 1),
        ]:
            g = sub[mask]
            if len(g) < 5:
                continue
            vals = g[WBC_PGS].to_numpy(float)
            mean, lo, hi = _boot_mean_ci(vals)
            rows.append(
                {
                    "super_pop": sp,
                    "duffy_null_gt": label,
                    "contrast": "dominant",
                    "n": len(g),
                    "mean_wbc_pgs": mean,
                    "sd_wbc_pgs": float(np.std(vals, ddof=1)) if len(g) > 1 else 0.0,
                    "mean_lo": lo,
                    "mean_hi": hi,
                    "underpowered": len(g) < 10,
                    "allele_definition": f"copies of {DUFFY_NULL_ALLELE} (Duffy-null)",
                }
            )
        # Primary AFR-style homozygous contrast: dose==2 vs dose<2
        for label, mask in [
            ("null_hom_2", sub["duffy_null_gt"] == 2),
            ("null_lt2", sub["duffy_null_gt"] < 2),
        ]:
            g = sub[mask]
            if len(g) < 3:
                continue
            vals = g[WBC_PGS].to_numpy(float)
            mean, lo, hi = _boot_mean_ci(vals)
            rows.append(
                {
                    "super_pop": sp,
                    "duffy_null_gt": label,
                    "contrast": "homozygous_vs_lt2",
                    "n": len(g),
                    "mean_wbc_pgs": mean,
                    "sd_wbc_pgs": float(np.std(vals, ddof=1)) if len(g) > 1 else 0.0,
                    "mean_lo": lo,
                    "mean_hi": hi,
                    "underpowered": len(g) < 10,
                    "allele_definition": f"copies of {DUFFY_NULL_ALLELE} (Duffy-null)",
                }
            )
        g2 = sub[sub["duffy_null_gt"] == 2]
        g0 = sub[sub["duffy_null_gt"] < 2]
        if len(g2) >= 5 and len(g0) >= 3:
            d = float(g2[WBC_PGS].mean() - g0[WBC_PGS].mean())
            rng = np.random.default_rng(719)
            boots = []
            a = g2[WBC_PGS].to_numpy(float)
            b = g0[WBC_PGS].to_numpy(float)
            for _ in range(400):
                boots.append(
                    float(np.mean(a[rng.integers(0, len(a), len(a))]) - np.mean(b[rng.integers(0, len(b), len(b))]))
                )
            rows.append(
                {
                    "super_pop": sp,
                    "duffy_null_gt": "delta_hom2_minus_lt2",
                    "contrast": "homozygous_vs_lt2",
                    "n": len(g2) + len(g0),
                    "mean_wbc_pgs": d,
                    "sd_wbc_pgs": np.nan,
                    "mean_lo": float(np.percentile(boots, 2.5)),
                    "mean_hi": float(np.percentile(boots, 97.5)),
                    "underpowered": len(g0) < 10,
                    "allele_definition": f"copies of {DUFFY_NULL_ALLELE} (Duffy-null)",
                    "n_hom2": len(g2),
                    "n_lt2": len(g0),
                }
            )
    return pd.DataFrame(rows)


def af_by_superpop(dose: pd.DataFrame, panel: pd.DataFrame) -> pd.DataFrame:
    m = dose.merge(panel[["sample_id", "super_pop"]], on="sample_id", how="inner")
    m = m.dropna(subset=["duffy_null_dose", "super_pop"])
    rows = []
    for sp, g in m.groupby("super_pop"):
        af = float(g["duffy_null_dose"].mean() / 2.0)
        rows.append(
            {
                "super_pop": sp,
                "n": len(g),
                "duffy_null_af": af,
                "expected_note": "AFR high (~0.8+); EUR near 0 for Duffy-null C",
            }
        )
    return pd.DataFrame(rows)


def _variant_pos(vid: str) -> tuple[str, int] | None:
    parts = str(vid).replace("chr", "").split(":")
    if len(parts) < 2:
        return None
    try:
        return parts[0], int(parts[1])
    except Exception:
        return None


def partition_weights(weights: pd.DataFrame) -> dict[str, pd.DataFrame]:
    rows_snp, rows_win, rows_out = [], [], []
    for _, r in weights.iterrows():
        vp = _variant_pos(r["variant_id"])
        if vp is None:
            rows_out.append(r)
            continue
        chrom, pos = vp
        if chrom == str(DUFFY_CHR) and pos == DUFFY_POS:
            rows_snp.append(r)
            rows_win.append(r)
        elif chrom == str(DUFFY_CHR) and abs(pos - DUFFY_POS) <= DUFFY_WINDOW_BP:
            rows_win.append(r)
        else:
            rows_out.append(r)
    return {
        "snp": pd.DataFrame(rows_snp),
        "window": pd.DataFrame(rows_win),
        "outside": pd.DataFrame(rows_out),
    }


def plink_score_weights(weights: pd.DataFrame, work: Path, tag: str) -> pd.DataFrame:
    """Score weights on available autosomal pfiles; sum per sample."""
    if weights.empty:
        return pd.DataFrame(columns=["sample_id", f"score_{tag}"])
    work.mkdir(parents=True, exist_ok=True)
    wpath = work / f"weights_{tag}.tsv"
    # PLINK2 --score: ID, allele, weight
    out_w = weights.copy()
    out_w["variant_id"] = out_w["variant_id"].astype(str)
    out_w[["variant_id", "effect_allele", "effect_weight"]].to_csv(wpath, sep="\t", index=False, header=False)
    parts = []
    for chrom in range(1, 23):
        pfile = PFILE_ROOT / f"chr{chrom}.score"
        if not (Path(str(pfile) + ".pgen").exists() or Path(str(pfile) + ".bed").exists()):
            continue
        # Filter weights to this chrom
        sub = out_w[out_w["variant_id"].astype(str).str.replace("^chr", "", regex=True).str.startswith(f"{chrom}:")]
        if sub.empty:
            continue
        wchrom = work / f"weights_{tag}_chr{chrom}.tsv"
        sub[["variant_id", "effect_allele", "effect_weight"]].to_csv(wchrom, sep="\t", index=False, header=False)
        out = work / f"score_{tag}_chr{chrom}"
        cmd = [
            plink2_bin(),
            "--pfile",
            str(pfile),
            "--score",
            str(wchrom),
            "1",
            "2",
            "3",
            "header-read" if False else "no-mean-imputation",
            "--out",
            str(out),
        ]
        # simpler: cols without header
        cmd = [
            plink2_bin(),
            "--pfile",
            str(pfile),
            "--score",
            str(wchrom),
            "1",
            "2",
            "3",
            "no-mean-imputation",
            "--out",
            str(out),
        ]
        r = subprocess.run(cmd, capture_output=True, text=True)
        sscore = Path(str(out) + ".sscore")
        if r.returncode != 0 or not sscore.exists():
            continue
        sc = pd.read_csv(sscore, sep=r"\s+")
        id_col = "#IID" if "#IID" in sc.columns else ("IID" if "IID" in sc.columns else sc.columns[0])
        score_col = next((c for c in sc.columns if "SCORE" in c.upper() or c.endswith("_SUM")), sc.columns[-1])
        parts.append(pd.DataFrame({"sample_id": sc[id_col].astype(str), "part": pd.to_numeric(sc[score_col], errors="coerce")}))
    if not parts:
        return pd.DataFrame(columns=["sample_id", f"score_{tag}"])
    allp = pd.concat(parts, ignore_index=True)
    summed = allp.groupby("sample_id", as_index=False)["part"].sum()
    return summed.rename(columns={"part": f"score_{tag}"})


def ackr1_decomposition(scores: pd.DataFrame, dose: pd.DataFrame, work: Path) -> pd.DataFrame:
    weights = load_genomewide_weights(WBC_PGS)
    parts = partition_weights(weights)
    print(
        f"ACKR1 weight partition: snp={len(parts['snp'])} window={len(parts['window'])} outside={len(parts['outside'])}",
        flush=True,
    )
    snp_sc = plink_score_weights(parts["snp"], work, "snp")
    win_sc = plink_score_weights(parts["window"], work, "window")
    # Prefer total from score matrix; outside = total - window
    m = scores[["sample_id", "super_pop", WBC_PGS]].merge(dose, on="sample_id", how="inner")
    m = m.rename(columns={WBC_PGS: "score_total"})
    if not snp_sc.empty:
        m = m.merge(snp_sc, on="sample_id", how="left")
    else:
        m["score_snp"] = np.nan
    if not win_sc.empty:
        m = m.merge(win_sc, on="sample_id", how="left")
    else:
        m["score_window"] = np.nan
    m["score_snp"] = m.get("score_snp", np.nan)
    m["score_window"] = m.get("score_window", np.nan)
    m["score_outside"] = m["score_total"] - m["score_window"]
    m["score_window_rest"] = m["score_window"] - m["score_snp"]
    m["duffy_null_gt"] = m["duffy_null_dose"].round().clip(0, 2).astype(int)

    rows = []
    for sp in ["AFR", "EUR"]:
        sub = m[m["super_pop"] == sp].dropna(subset=["score_total"])
        g2 = sub[sub["duffy_null_gt"] == 2]
        g0 = sub[sub["duffy_null_gt"] < 2]
        if len(g2) < 5 or len(g0) < 3:
            rows.append(
                {
                    "super_pop": sp,
                    "status": "underpowered_contrast",
                    "n_hom2": len(g2),
                    "n_lt2": len(g0),
                }
            )
            continue
        rng = np.random.default_rng(719)

        def delta(col: str) -> tuple[float, float, float]:
            a = g2[col].to_numpy(float)
            b = g0[col].to_numpy(float)
            a = a[np.isfinite(a)]
            b = b[np.isfinite(b)]
            if len(a) < 3 or len(b) < 3:
                return float("nan"), float("nan"), float("nan")
            d0 = float(np.mean(a) - np.mean(b))
            boots = []
            for _ in range(400):
                boots.append(
                    float(np.mean(a[rng.integers(0, len(a), len(a))]) - np.mean(b[rng.integers(0, len(b), len(b))]))
                )
            return d0, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))

        d_tot, lo_t, hi_t = delta("score_total")
        d_win, lo_w, hi_w = delta("score_window")
        d_snp, lo_s, hi_s = delta("score_snp")
        d_out, lo_o, hi_o = delta("score_outside")
        d_wr, lo_wr, hi_wr = delta("score_window_rest")
        frac = d_win / d_tot if np.isfinite(d_win) and np.isfinite(d_tot) and abs(d_tot) > 1e-12 else np.nan
        # bootstrap fraction
        a_t = g2["score_total"].to_numpy(float)
        b_t = g0["score_total"].to_numpy(float)
        a_w = g2["score_window"].to_numpy(float)
        b_w = g0["score_window"].to_numpy(float)
        fracs = []
        for _ in range(400):
            ia = rng.integers(0, len(a_t), len(a_t))
            ib = rng.integers(0, len(b_t), len(b_t))
            dt = float(np.nanmean(a_t[ia]) - np.nanmean(b_t[ib]))
            dw = float(np.nanmean(a_w[ia]) - np.nanmean(b_w[ib]))
            if abs(dt) > 1e-12:
                fracs.append(dw / dt)
        rows.append(
            {
                "super_pop": sp,
                "status": "ok",
                "contrast": "hom2_minus_lt2",
                "n_hom2": len(g2),
                "n_lt2": len(g0),
                "delta_total": d_tot,
                "delta_total_lo": lo_t,
                "delta_total_hi": hi_t,
                "delta_window": d_win,
                "delta_window_lo": lo_w,
                "delta_window_hi": hi_w,
                "delta_snp": d_snp,
                "delta_snp_lo": lo_s,
                "delta_snp_hi": hi_s,
                "delta_window_rest": d_wr,
                "delta_outside": d_out,
                "delta_outside_lo": lo_o,
                "delta_outside_hi": hi_o,
                "ackr1_fraction": frac,
                "ackr1_fraction_lo": float(np.percentile(fracs, 2.5)) if fracs else np.nan,
                "ackr1_fraction_hi": float(np.percentile(fracs, 97.5)) if fracs else np.nan,
                "n_snp_weights": len(parts["snp"]),
                "n_window_weights": len(parts["window"]),
                "n_outside_weights": len(parts["outside"]),
            }
        )
    return pd.DataFrame(rows)


def risk_near_duffy(preds: Path) -> dict:
    if not preds.exists():
        return {"status": "missing_predictions"}
    try:
        import pyarrow.parquet as pq

        schema = pq.read_schema(preds)
        cols = [c for c in ["variant_id", "predicted_risk", "risk", "y_pred"] if c in schema.names]
        df = pd.read_parquet(preds, columns=cols)
    except Exception:
        df = pd.read_parquet(preds)
    if "predicted_risk" not in df.columns:
        alt = next((c for c in ["risk", "y_pred"] if c in df.columns), None)
        if alt is None:
            return {"status": "bad_columns"}
        df = df.rename(columns={alt: "predicted_risk"})
    df = df.groupby("variant_id", as_index=False)["predicted_risk"].mean()
    parts = df["variant_id"].astype(str).str.split(":", n=3, expand=True)
    df["chrom"] = parts[0].str.replace("^chr", "", regex=True)
    df["pos"] = pd.to_numeric(parts[1], errors="coerce")
    chr1 = df[df["chrom"] == str(DUFFY_CHR)].dropna(subset=["pos"])
    near = chr1[(chr1["pos"] - DUFFY_POS).abs() <= DUFFY_WINDOW_BP]
    far = chr1[(chr1["pos"] - DUFFY_POS).abs() > DUFFY_WINDOW_BP]
    return {
        "status": "ok",
        "n_near": int(len(near)),
        "n_far_chr1": int(len(far)),
        "mean_risk_near": float(near["predicted_risk"].mean()) if len(near) else float("nan"),
        "mean_risk_far_chr1": float(far["predicted_risk"].mean()) if len(far) else float("nan"),
        "note": "Portability-risk near Duffy is not the positive control; genotype–score is.",
    }


def plot_strata(strata: pd.DataFrame, fig_dir: Path, audit: dict) -> Path | None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None
    fig_dir.mkdir(parents=True, exist_ok=True)
    num = strata[pd.to_numeric(strata["duffy_null_gt"], errors="coerce").notna()].copy()
    num["duffy_null_gt"] = num["duffy_null_gt"].astype(int)
    if "underpowered" in num.columns:
        num = num[~num["underpowered"].astype(bool)]
    fig, ax = plt.subplots(figsize=(7, 4))
    for sp, color in [("AFR", "#c45c26"), ("EUR", "#2c5f7c")]:
        sub = num[num["super_pop"] == sp]
        if sub.empty:
            continue
        ax.errorbar(
            sub["duffy_null_gt"],
            sub["mean_wbc_pgs"],
            yerr=[sub["mean_wbc_pgs"] - sub["mean_lo"], sub["mean_hi"] - sub["mean_wbc_pgs"]],
            marker="o",
            label=sp,
            color=color,
            capsize=3,
        )
    ax.set_xlabel(f"Copies of Duffy-null allele ({audit.get('duffy_null_allele_target', 'C')})")
    ax.set_ylabel(f"Mean {WBC_PGS} score")
    ax.set_title(
        f"{DUFFY_RSID} {audit.get('ref','?')}/{audit.get('alt','?')} "
        f"(counted {audit.get('plink_counted_allele','?')}; "
        f"flipped={audit.get('dosage_flipped_to_null')})"
    )
    ax.legend()
    ax.set_xticks([0, 1, 2])
    out = fig_dir / "fig_duffy_wbc_genotype.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    work = ROOT / "data/interim/duffy_control"
    print("Extracting Duffy dosages with allele audit ...", flush=True)
    dose, audit = extract_duffy_dosages(Path(args.pfile), work)
    (out_dir / "duffy_allele_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2), flush=True)

    scores = pd.read_parquet(args.score_matrix)
    if WBC_PGS not in scores.columns:
        raise SystemExit(f"{WBC_PGS} missing from score matrix")
    if "super_pop" not in scores.columns and PANEL.exists():
        panel = pd.read_parquet(PANEL)
        scores = scores.merge(panel[["sample_id", "super_pop"]], on="sample_id", how="left")

    panel = pd.read_parquet(PANEL) if PANEL.exists() else scores[["sample_id", "super_pop"]]
    af = af_by_superpop(dose, panel)
    af.to_csv(out_dir / "duffy_null_af_by_superpop.csv", index=False, float_format="%.6f")

    strata = genotype_strata(scores, dose)
    strata.to_csv(out_dir / "duffy_wbc_genotype_strata.csv", index=False, float_format="%.6f")
    pd.DataFrame([audit]).to_csv(out_dir / "duffy_allele_audit.csv", index=False)

    print("ACKR1 score decomposition ...", flush=True)
    decomp = ackr1_decomposition(scores, dose, work / "ackr1_scores")
    decomp.to_csv(out_dir / "duffy_ackr1_score_decomposition.csv", index=False, float_format="%.6g")
    print(decomp.to_string(index=False), flush=True)

    fig = plot_strata(strata, Path(args.fig_dir), audit)
    if fig:
        print(f"Saved {fig}", flush=True)

    gate = risk_near_duffy(Path(args.predictions))
    pd.DataFrame([gate]).to_csv(out_dir / "duffy_positive_control_genomewide.csv", index=False)
    print("DUFFY_CONTROL_DONE", flush=True)


if __name__ == "__main__":
    main()
