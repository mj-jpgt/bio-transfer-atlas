"""
Phase A4: Subpopulation-resolved AF features from existing pop_*.keep files.

Computes within-superpop AF variance and max pairwise subpop AF diff per variant.
Fixes plink2 .afreq parsing (#CHROM header must not be treated as a comment line).
"""
from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from bta_plink import plink2_bin  # noqa: E402

PLINK2 = plink2_bin()
INTERIM = ROOT / "data/interim/1000g_grch38"
KEEP = ROOT / "data/features/af/keep_files"
OUT_DIR = ROOT / "data/features/af"

SUPERPOP_TO_POPS = {
    "AFR": ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"],
    "EUR": ["CEU", "TSI", "FIN", "GBR", "IBS"],
    "EAS": ["CHB", "JPT", "CHS", "CDX", "KHV"],
    "SAS": ["GIH", "PJL", "BEB", "STU", "ITU"],
    "AMR": ["MXL", "PUR", "CLM", "PEL"],
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chroms", default="22", help="Comma list e.g. 1,2,22 or 1-22")
    p.add_argument("--memory-mb", type=int, default=int(os.environ.get("BTA_PLINK_MEMORY_MB", "4096")))
    p.add_argument("--jobs", type=int, default=1, help="Parallel chromosomes")
    return p.parse_args()


def parse_chroms(spec: str) -> list[str]:
    s = spec.strip()
    if "-" in s and "," not in s:
        a, b = s.split("-", 1)
        return [str(x) for x in range(int(a), int(b) + 1)]
    return [c.strip() for c in s.split(",") if c.strip()]


def read_afreq(path: Path) -> pd.DataFrame:
    """Parse plink2 .afreq; header may be #CHROM ... ALT_FREQS ..."""
    rows = []
    header = None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                # plink2 header: #CHROM ID REF ALT ALT_FREQS OBS_CT
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if header is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            rows.append(dict(zip(header, parts)))
    if not rows:
        return pd.DataFrame(columns=["variant_id", "af"])
    df = pd.DataFrame(rows)
    id_col = "ID" if "ID" in df.columns else df.columns[1]
    freq_col = "ALT_FREQS" if "ALT_FREQS" in df.columns else None
    if freq_col is None:
        # never fall back to OBS_CT
        for c in df.columns:
            if "FREQ" in c.upper():
                freq_col = c
                break
    if freq_col is None:
        raise RuntimeError(f"No ALT_FREQS in {path}; cols={list(df.columns)}")
    out = df[[id_col, freq_col]].rename(columns={id_col: "variant_id", freq_col: "af"})
    out["af"] = pd.to_numeric(out["af"], errors="coerce")
    # allele frequencies must be in [0, 1]
    out = out[(out["af"] >= 0) & (out["af"] <= 1)]
    return out


def plink_freq(pfile: Path, keep: Path, out_prefix: Path, memory_mb: int) -> pd.DataFrame | None:
    import subprocess

    if not keep.exists() or keep.stat().st_size == 0:
        return None
    cmd = [
        PLINK2,
        "--pfile",
        str(pfile),
        "--keep",
        str(keep),
        "--freq",
        "--memory",
        str(memory_mb),
        "--out",
        str(out_prefix),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    afreq = Path(str(out_prefix) + ".afreq")
    if r.returncode != 0 or not afreq.exists():
        return None
    return read_afreq(afreq)


def process_chrom(chrom: str, memory_mb: int) -> Path | None:
    pfile = INTERIM / f"chr{chrom}.qc"
    if not (INTERIM / f"chr{chrom}.qc.pgen").exists():
        pfile = INTERIM / f"chr{chrom}.score"
    if not Path(str(pfile) + ".pgen").exists():
        print(f"  chr{chrom}: no pgen, skip", flush=True)
        return None

    work = OUT_DIR / "subpop_tmp" / f"chr{chrom}"
    work.mkdir(parents=True, exist_ok=True)
    frames = {}
    for sp, pops in SUPERPOP_TO_POPS.items():
        afs = []
        for pop in pops:
            keep = KEEP / f"pop_{pop}.keep"
            pref = work / f"{pop}"
            df = plink_freq(pfile, keep, pref, memory_mb)
            if df is not None and len(df):
                afs.append(df.rename(columns={"af": pop}).set_index("variant_id"))
        if len(afs) < 2:
            continue
        mat = pd.concat(afs, axis=1)
        # drop rows with <2 non-null pops
        valid = mat.dropna(axis=0, thresh=2)
        if valid.empty:
            continue
        frames[f"AF_var_within_{sp}"] = valid.var(axis=1, skipna=True)
        frames[f"AF_maxdiff_within_{sp}"] = valid.max(axis=1) - valid.min(axis=1)

    if not frames:
        print(f"  chr{chrom}: no subpop freqs", flush=True)
        return None
    out = pd.DataFrame(frames)
    out.index.name = "variant_id"
    out = out.reset_index()
    within_cols = [c for c in out.columns if c.startswith("AF_var_within_")]
    out["AF_var_within_mean"] = out[within_cols].mean(axis=1)
    out["AF_var_within_max"] = out[within_cols].max(axis=1)

    # sanity: maxdiff must be in [0, 1]
    for c in out.columns:
        if c.startswith("AF_maxdiff_"):
            bad = ((out[c] < 0) | (out[c] > 1)).sum()
            if bad:
                raise RuntimeError(f"chr{chrom}: {bad} bad values in {c} (AF parse broken)")

    path = OUT_DIR / f"subpop_af_features.chr{chrom}.parquet"
    out.to_parquet(path, index=False)
    print(
        f"  chr{chrom}: saved {path} ({len(out):,} rows); "
        f"AF_maxdiff_within_AFR mean={out.get('AF_maxdiff_within_AFR', pd.Series([np.nan])).mean():.4f}",
        flush=True,
    )
    return path


def main() -> None:
    args = parse_args()
    chroms = parse_chroms(args.chroms)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not KEEP.exists():
        raise SystemExit(f"Missing keep files at {KEEP}")

    paths: list[Path] = []
    jobs = max(1, args.jobs)
    if jobs == 1:
        for chrom in chroms:
            print(f"=== chr{chrom} ===", flush=True)
            p = process_chrom(chrom, args.memory_mb)
            if p:
                paths.append(p)
    else:
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(process_chrom, c, args.memory_mb): c for c in chroms}
            for fut in as_completed(futs):
                chrom = futs[fut]
                print(f"=== chr{chrom} done ===", flush=True)
                p = fut.result()
                if p:
                    paths.append(p)

    status_rows = []
    for chrom in chroms:
        path = OUT_DIR / f"subpop_af_features.chr{chrom}.parquet"
        if path.exists():
            status_rows.append({"chrom": chrom, "status": "ok", "path": str(path)})
        else:
            status_rows.append({"chrom": chrom, "status": "failed_or_skipped", "path": ""})
    st = ROOT / "results/tables/subpop_af_features_status.csv"
    st.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(status_rows).to_csv(st, index=False)
    print(f"Saved {st}")

    if not paths:
        raise SystemExit("No subpop features produced")

    rows = []
    for p in paths:
        chrom = p.name.split("chr")[1].split(".")[0]
        sub = pd.read_parquet(p)
        af_path = OUT_DIR / f"af_divergence_features.chr{chrom}.parquet"
        if af_path.exists():
            af = pd.read_parquet(af_path, columns=["variant_id", "AF_var"])
            m = sub.merge(af, on="variant_id", how="inner")
            if len(m) and "AF_var" in m.columns:
                rows.append(
                    {
                        "chrom": chrom,
                        "mean_AF_var_between_superpop": float(m["AF_var"].mean()),
                        "mean_AF_var_within_superpop": float(m["AF_var_within_mean"].mean()),
                        "ratio_within_over_between": float(
                            m["AF_var_within_mean"].mean() / max(m["AF_var"].mean(), 1e-12)
                        ),
                        "n": len(m),
                    }
                )
    if rows:
        summary = pd.DataFrame(rows)
        spath = ROOT / "results/tables/subpop_af_variance_decomposition.csv"
        summary.to_csv(spath, index=False, float_format="%.6g")
        print(summary.to_string(index=False))
        print(f"Saved {spath}")


if __name__ == "__main__":
    main()
