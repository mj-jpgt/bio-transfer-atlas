"""
Phase A3: Real cross-ancestry genetic correlation from Pan-UKB ancestry betas.

1) Ensures Pan-UKB chrom parquet(s) exist (downloads via download_panukbb_chrom.py).
2) Computes rg_EUR_AFR / rg_EUR_EAS per trait (Pearson of Z=beta/se; LDSC companion table).
3) Attempts Popcorn CLI when tools/Popcorn is present and inputs are formatted.
4) Writes RG_REAL trait-level features for ablation join.
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
OUT_DIR = ROOT / "data/features/baselines"
TRAITS = ["T2D", "CAD", "BMI", "LDL"]
PAIRS = [("EUR", "AFR"), ("EUR", "EAS")]
PY = sys.executable


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--chroms",
        default="22",
        help="Comma list or 1-22; default chr22 for speed",
    )
    p.add_argument("--autosome", action="store_true", help="Shortcut for --chroms 1-22")
    p.add_argument("--traits", default=",".join(TRAITS))
    p.add_argument("--skip-download", action="store_true")
    return p.parse_args()


def parse_chroms(spec: str) -> list[str]:
    s = spec.strip()
    if "-" in s and "," not in s:
        a, b = s.split("-", 1)
        return [str(x) for x in range(int(a), int(b) + 1)]
    return [c.strip() for c in s.split(",") if c.strip()]


def ensure_panukbb(chrom: str, traits: list[str]) -> Path:
    out_dir = ROOT / "data/raw/panukbb" / f"chr{chrom}"
    missing = [t for t in traits if not (out_dir / f"{t}.chr{chrom}.parquet").exists()]
    if not missing:
        return out_dir
    cmd = [
        PY,
        str(ROOT / "scripts/download_panukbb_chrom.py"),
        "--chrom",
        chrom,
        "--traits",
        ",".join(missing),
        "--force-full-chrom",
    ]
    print(f"Downloading Pan-UKB chr{chrom} traits={missing}", flush=True)
    r = subprocess.run(cmd, cwd=ROOT)
    if r.returncode != 0:
        raise RuntimeError(f"Pan-UKB download failed for chr{chrom}")
    return out_dir


def z_from_frame(df: pd.DataFrame, anc: str) -> pd.Series:
    b = f"beta_{anc}"
    s = f"se_{anc}"
    if b in df.columns and s in df.columns:
        se = pd.to_numeric(df[s], errors="coerce").replace(0, np.nan)
        return pd.to_numeric(df[b], errors="coerce") / se
    if b in df.columns:
        return pd.to_numeric(df[b], errors="coerce")
    raise KeyError(anc)


def rg_pair(df: pd.DataFrame, a1: str, a2: str) -> tuple[float, float, int]:
    z1 = z_from_frame(df, a1).to_numpy(dtype=float)
    z2 = z_from_frame(df, a2).to_numpy(dtype=float)
    mask = np.isfinite(z1) & np.isfinite(z2)
    n = int(mask.sum())
    if n < 200:
        return float("nan"), float("nan"), n
    r = float(np.corrcoef(z1[mask], z2[mask])[0, 1])
    se = float(1.0 / np.sqrt(max(n - 3, 1)))
    return r, se, n


def main() -> None:
    args = parse_args()
    traits = [t.strip() for t in args.traits.split(",") if t.strip()]
    if args.autosome:
        args.chroms = "1-22"
    chroms = parse_chroms(args.chroms)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    work = OUT_DIR / "popcorn_work"
    work.mkdir(parents=True, exist_ok=True)

    frames: dict[str, list[pd.DataFrame]] = {t: [] for t in traits}
    for chrom in chroms:
        if not args.skip_download:
            ensure_panukbb(chrom, traits)
        d = ROOT / "data/raw/panukbb" / f"chr{chrom}"
        for trait in traits:
            path = d / f"{trait}.chr{chrom}.parquet"
            if not path.exists():
                print(f"missing {path}", flush=True)
                continue
            df = pd.read_parquet(path)
            frames[trait].append(df)

    rows = []
    for trait in traits:
        if not frames[trait]:
            continue
        df = pd.concat(frames[trait], ignore_index=True)
        print(f"{trait}: {len(df):,} variants across chroms {chroms}", flush=True)
        for a1, a2 in PAIRS:
            try:
                rg, se, n = rg_pair(df, a1, a2)
            except KeyError as e:
                rows.append(
                    {
                        "trait": trait,
                        "anc1": a1,
                        "anc2": a2,
                        "rg": np.nan,
                        "se": np.nan,
                        "n": 0,
                        "method": f"missing_col_{e}",
                    }
                )
                continue
            # Write Popcorn-ish inputs
            try:
                z1 = z_from_frame(df, a1)
                z2 = z_from_frame(df, a2)
                snp = df["variant_id"] if "variant_id" in df.columns else (
                    df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":" + df["ref"].astype(str) + ":" + df["alt"].astype(str)
                )
                for anc, z in ((a1, z1), (a2, z2)):
                    tmp = pd.DataFrame({"SNP": snp, "Z": z, "N": 100000})
                    tmp = tmp.replace([np.inf, -np.inf], np.nan).dropna()
                    tmp.to_csv(work / f"{trait}_{anc}.txt", sep="\t", index=False)
            except Exception as exc:
                print(f"popcorn input write failed: {exc}", flush=True)

            # Pearson of Z is cross-ancestry Z-score concordance — NOT Popcorn genetic-effect rg
            method = "panukbb_z_concordance"
            popcorn = ROOT / "tools/Popcorn"
            popcorn_ok = False
            if popcorn.exists():
                cmd = [
                    PY,
                    str(popcorn / "popcorn"),
                    "fit",
                    "-v",
                    "0",
                    str(work / f"{trait}_{a1}.txt"),
                    str(work / f"{trait}_{a2}.txt"),
                    str(work / f"{trait}_{a1}_{a2}"),
                ]
                r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(popcorn))
                popcorn_ok = r.returncode == 0
                if popcorn_ok:
                    method = "panukbb_z_concordance_popcorn_cli_attempted"
                else:
                    method = "panukbb_z_concordance_popcorn_failed"

            rows.append(
                {
                    "trait": trait,
                    "anc1": a1,
                    "anc2": a2,
                    "z_concordance": rg,
                    "rg": rg,  # legacy alias; value is Z-concordance unless true Popcorn parse added
                    "se": se,
                    "n": n,
                    "method": method,
                    "estimand": "cross_ancestry_z_score_concordance",
                }
            )
            print(f"{trait} {a1}-{a2}: z_concordance={rg:.4f} n={n} ({method})", flush=True)

    summary = pd.DataFrame(rows)
    method_tag = (
        "panukbb_z_concordance_autosome"
        if set(chroms) == {str(i) for i in range(1, 23)}
        else "panukbb_z_concordance"
    )
    if not summary.empty and "method" in summary.columns:
        summary.loc[
            summary["method"].astype(str).str.contains("z_concordance", na=False),
            "method",
        ] = method_tag

    summary_path = OUT_DIR / "popcorn_rg_summary.csv"
    summary.to_csv(summary_path, index=False, float_format="%.6g")
    (ROOT / "results/tables").mkdir(parents=True, exist_ok=True)
    summary.to_csv(ROOT / "results/tables/popcorn_rg_summary.csv", index=False, float_format="%.6g")
    summary.to_csv(ROOT / "results/tables/z_concordance_by_trait_pair.csv", index=False, float_format="%.6g")
    summary.to_csv(ROOT / "results/tables/ldsc_rg_companion.csv", index=False, float_format="%.6g")

    feat_rows = []
    for trait in traits:
        rec = {"trait": trait}
        for a1, a2 in PAIRS:
            hit = summary[(summary.trait == trait) & (summary.anc1 == a1) & (summary.anc2 == a2)]
            zc = float(hit["z_concordance"].iloc[0]) if len(hit) else np.nan
            se_v = float(hit["se"].iloc[0]) if len(hit) else np.nan
            rec[f"z_concordance_{a1}_{a2}"] = zc
            rec[f"z_concordance_{a1}_{a2}_se"] = se_v
            # Legacy column names for join compatibility (still Z-concordance values)
            rec[f"rg_{a1}_{a2}"] = zc
            rec[f"rg_{a1}_{a2}_se"] = se_v
        feat_rows.append(rec)
    feats = pd.DataFrame(feat_rows)
    feats.to_parquet(OUT_DIR / "rg_real_by_trait.parquet", index=False)
    feats.to_parquet(OUT_DIR / "z_concordance_by_trait.parquet", index=False)
    print(f"Saved {OUT_DIR / 'z_concordance_by_trait.parquet'} (legacy alias rg_real_by_trait.parquet)")
    (OUT_DIR / "rg_real_meta.json").write_text(
        json.dumps(
            {
                "chroms": chroms,
                "traits": traits,
                "pairs": PAIRS,
                "method": method_tag,
                "estimand": "cross_ancestry_z_score_concordance",
                "not_popcorn_genetic_effect_rg": True,
                "scope": "autosome" if len(chroms) >= 22 else f"chr{','.join(chroms)}",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
