"""
Compute LD-divergence features for a target chromosome.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PLINK2 = str(ROOT / "tools/plink2/plink2.exe")
DEFAULT_PLINK_MEMORY_MB = 2048

INTERIM = ROOT / "data/interim/1000g_grch38"
AF_OUT = ROOT / "data/features/af"
LD_OUT = ROOT / "data/features/ld"
LD_OUT.mkdir(parents=True, exist_ok=True)

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
KEEP_DIR = AF_OUT / "keep_files"
LD_WINDOW_KB = 500
THREADS = 8


def normalize_chr(value: str) -> str:
    s = str(value).strip()
    return s[3:] if s.lower().startswith("chr") else s


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22", help="Chromosome number/name (e.g. 22 or chr22)")
    p.add_argument("--threads", type=int, default=THREADS)
    p.add_argument("--memory-mb", type=int, default=DEFAULT_PLINK_MEMORY_MB)
    return p.parse_args()


def run(cmd: list[str], desc: str) -> None:
    print(f"  [{desc}]", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        for line in r.stderr.strip().split("\n")[-10:]:
            if line.strip():
                print(f"    {line}")
        raise RuntimeError(f"FAILED: {desc}")


def maybe_write_legacy(chrom: str, src: Path, legacy_name: str) -> None:
    if chrom == "22":
        pd.read_parquet(src).to_parquet(LD_OUT / legacy_name, index=False)


def compute_ld_scores(
    sp: str, chrom: str, extract_file: Path, threads: int, memory_mb: int
) -> pd.DataFrame:
    suffix = f".chr{chrom}"
    cache_path = LD_OUT / f"ld_scores_{sp}{suffix}.parquet"
    if cache_path.exists():
        print(f"  {sp}: LD scores cached, loading ...")
        return pd.read_parquet(cache_path)

    vcor_prefix = LD_OUT / f"r2_{sp}{suffix}"
    vcor_file = LD_OUT / f"r2_{sp}{suffix}.vcor"
    if not vcor_file.exists():
        run(
            [
                PLINK2,
                "--memory",
                str(memory_mb),
                "--pfile",
                str(INTERIM / f"chr{chrom}.score"),
                "--extract",
                str(extract_file),
                "--keep",
                str(KEEP_DIR / f"{sp}.keep"),
                "--r2-unphased",
                "yes-really",
                "cols=id",
                "--ld-window-kb",
                str(LD_WINDOW_KB),
                "--ld-window-r2",
                "0.01",
                "--threads",
                str(threads),
                "--out",
                str(vcor_prefix),
            ],
            f"r2-unphased chr{chrom} {sp}",
        )

    print(f"  {sp}: reading r² table in chunks ...", flush=True)
    sums: dict[str, float] = {}
    chunk_iter = pd.read_csv(
        vcor_file,
        sep="\t",
        usecols=["#ID_A", "ID_B", "UNPHASED_R2"],
        dtype={"#ID_A": str, "ID_B": str, "UNPHASED_R2": np.float32},
        chunksize=2_000_000,
    )
    for i, chunk in enumerate(chunk_iter):
        chunk = chunk.rename(columns={"#ID_A": "v1", "ID_B": "v2", "UNPHASED_R2": "R2"})
        for row in chunk.itertuples(index=False):
            sums[row.v1] = sums.get(row.v1, 0.0) + row.R2
            sums[row.v2] = sums.get(row.v2, 0.0) + row.R2
        if i % 20 == 0:
            print(f"    {sp}: {(i + 1) * 2_000_000:,} rows processed ...", flush=True)
    ld = pd.DataFrame(list(sums.items()), columns=["variant_id", f"LD_{sp}"])
    ld.to_parquet(cache_path, index=False)
    print(f"  {sp}: {len(ld):,} variants  median LD = {ld[f'LD_{sp}'].median():.2f}")
    return ld


def main() -> None:
    args = parse_args()
    chrom = normalize_chr(args.chrom)
    suffix = f".chr{chrom}"
    extract_file = AF_OUT / f"pgs_variants{suffix}.extract"
    if not extract_file.exists() and chrom == "22":
        extract_file = AF_OUT / "pgs_variants.extract"
    if not extract_file.exists():
        raise FileNotFoundError(
            f"Missing {extract_file}. Run compute_af_features.py --chrom {chrom} first."
        )
    print(f"Matched PGS variants (targets) chr{chrom}: {sum(1 for _ in open(extract_file, encoding='utf-8')):,}")

    sp_ld: dict[str, pd.DataFrame] = {}
    for sp in SUPERPOPS:
        sp_ld[sp] = compute_ld_scores(sp, chrom, extract_file, args.threads, args.memory_mb)

    print("\nMerging LD scores ...")
    all_variants = pd.read_csv(extract_file, header=None, names=["variant_id"], dtype=str)
    merged = all_variants.copy()
    for sp in SUPERPOPS:
        merged = merged.merge(sp_ld[sp], on="variant_id", how="left")

    ld_cols = [f"LD_{sp}" for sp in SUPERPOPS]
    merged[ld_cols] = merged[ld_cols].fillna(0.0)
    ld_mat = merged[ld_cols].values.astype(float)

    merged["LD_max_diff"] = np.nanmax(ld_mat, axis=1) - np.nanmin(ld_mat, axis=1)
    merged["LD_EUR_AFR_ratio"] = merged["LD_EUR"] / merged["LD_AFR"].replace(0, np.nan)
    merged["LD_EUR_EAS_ratio"] = merged["LD_EUR"] / merged["LD_EAS"].replace(0, np.nan)

    def ld_entropy(row: pd.Series) -> float:
        vals = np.array([row[c] for c in ld_cols], dtype=float)
        vals = np.clip(vals, 1e-9, None)
        vals = vals / vals.sum()
        return float(-np.sum(vals * np.log(vals)) / np.log(len(vals)))

    merged["LD_entropy"] = merged.apply(ld_entropy, axis=1)

    out = LD_OUT / f"ld_divergence_features{suffix}.parquet"
    merged.to_parquet(out, index=False)
    print(f"\nSaved: {out.name} ({len(merged):,} variants)")
    maybe_write_legacy(chrom, out, "ld_divergence_features.parquet")

    n_total = len(all_variants)
    n_with_ld = merged[["LD_AFR", "LD_EUR", "LD_EAS", "LD_SAS"]].notna().all(axis=1).sum()
    pct = n_with_ld / n_total * 100
    print(f"\nGate check chr{chrom}: {n_with_ld:,}/{n_total:,} variants with LD in AFR+EUR+EAS+SAS ({pct:.1f}%)")
    print("  PASS" if pct >= 70 else "  WARNING: <70%")


if __name__ == "__main__":
    main()
