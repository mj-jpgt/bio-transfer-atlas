"""
Compute AF-divergence features for a target chromosome.
"""
from __future__ import annotations

import argparse
import subprocess
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PLINK2 = str(ROOT / "tools/plink2/plink2.exe")
DEFAULT_PLINK_MEMORY_MB = 2048

INTERIM = ROOT / "data/interim/1000g_grch38"
META_DIR = ROOT / "data/raw/1000g/metadata"
PGS_OUT = ROOT / "data/processed/pgs_grch38"
AF_OUT = ROOT / "data/features/af"
AF_OUT.mkdir(parents=True, exist_ok=True)

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def normalize_chr(value: str) -> str:
    s = str(value).strip()
    return s[3:] if s.lower().startswith("chr") else s


def run(cmd: list[str], desc: str) -> None:
    print(f"  [{desc}]", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        for line in r.stderr.strip().split("\n")[-10:]:
            if line.strip():
                print(f"    {line}")
        raise RuntimeError(f"FAILED: {desc}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22", help="Chromosome number/name (e.g. 22 or chr22)")
    p.add_argument("--memory-mb", type=int, default=DEFAULT_PLINK_MEMORY_MB)
    return p.parse_args()


def collect_variant_ids(chrom: str) -> set[str]:
    print(f"Collecting matched variant IDs for chr{chrom} ...")
    all_variant_ids: set[str] = set()
    prefix = f"{chrom}:"
    for pgs_dir in PGS_OUT.iterdir():
        if not pgs_dir.is_dir():
            continue
        tsv = pgs_dir / f"{pgs_dir.name}.harmonized.tsv"
        if not tsv.exists():
            continue
        for chunk in pd.read_csv(tsv, sep="\t", usecols=["variant_id"], dtype=str, chunksize=250_000):
            vid = chunk["variant_id"].dropna()
            keep = vid[vid.str.startswith(prefix)]
            all_variant_ids.update(keep.tolist())
    print(f"  Total unique matched variants on chr{chrom}: {len(all_variant_ids):,}")
    return all_variant_ids


def write_keep_files() -> Path:
    panel_path = META_DIR / "integrated_call_samples_v3.20130502.ALL.panel"
    panel = pd.read_csv(panel_path, sep="\t")
    panel.columns = panel.columns.str.strip()
    rename = {}
    for c in panel.columns:
        lc = c.lower().strip()
        if lc in ("sample", "#sample", "individual_id"):
            rename[c] = "sample_id"
        elif lc == "pop":
            rename[c] = "pop"
        elif lc in ("super_pop", "super population code"):
            rename[c] = "super_pop"
    panel = panel.rename(columns=rename)

    keep_dir = AF_OUT / "keep_files"
    keep_dir.mkdir(exist_ok=True)
    for sp in SUPERPOPS:
        samp = panel[panel["super_pop"] == sp]["sample_id"].tolist()
        (keep_dir / f"{sp}.keep").write_text("\n".join(samp), encoding="utf-8")
    for pop in panel["pop"].dropna().unique():
        samp = panel[panel["pop"] == pop]["sample_id"].tolist()
        (keep_dir / f"pop_{pop}.keep").write_text("\n".join(samp), encoding="utf-8")
    return keep_dir


def maybe_write_legacy(chrom: str, src: Path, legacy_name: str) -> None:
    if chrom == "22":
        pd.read_parquet(src).to_parquet(AF_OUT / legacy_name, index=False)


def main() -> None:
    args = parse_args()
    chrom = normalize_chr(args.chrom)
    plink_mem = args.memory_mb
    suffix = f".chr{chrom}"
    pfile_prefix = INTERIM / f"chr{chrom}.score"
    if not Path(str(pfile_prefix) + ".pgen").exists():
        raise FileNotFoundError(
            f"Missing {pfile_prefix}.pgen. Run preprocess for chr{chrom} first."
        )

    all_variant_ids = collect_variant_ids(chrom)
    if not all_variant_ids:
        raise RuntimeError(f"No matched variants found for chr{chrom}")

    extract_file = AF_OUT / f"pgs_variants{suffix}.extract"
    extract_file.write_text("\n".join(sorted(all_variant_ids)), encoding="utf-8")
    if chrom == "22":
        (AF_OUT / "pgs_variants.extract").write_text(
            "\n".join(sorted(all_variant_ids)), encoding="utf-8"
        )

    keep_dir = write_keep_files()

    print("\nComputing AF per superpopulation ...")
    sp_afreqs: dict[str, pd.DataFrame] = {}
    for sp in SUPERPOPS:
        out_prefix = str(AF_OUT / f"freq_{sp}{suffix}")
        if Path(out_prefix + ".afreq").exists():
            print(f"  {sp}: freq file exists, loading ...")
        else:
            run(
                [
                    PLINK2,
                    "--memory",
                    str(plink_mem),
                    "--pfile",
                    str(pfile_prefix),
                    "--extract",
                    str(extract_file),
                    "--keep",
                    str(keep_dir / f"{sp}.keep"),
                    "--freq",
                    "--out",
                    out_prefix,
                ],
                f"freq chr{chrom} {sp}",
            )
        afreq = pd.read_csv(out_prefix + ".afreq", sep="\t")
        afreq = afreq.rename(
            columns={"#CHROM": "chr", "ID": "variant_id", "ALT_FREQS": f"AF_{sp}"}
        )
        sp_afreqs[sp] = afreq[["variant_id", f"AF_{sp}", "OBS_CT"]].copy()
        print(f"  {sp}: {len(afreq):,} variants")

    print("\nMerging superpop AFs ...")
    merged = sp_afreqs[SUPERPOPS[0]][["variant_id", f"AF_{SUPERPOPS[0]}"]].copy()
    for sp in SUPERPOPS[1:]:
        merged = merged.merge(sp_afreqs[sp][["variant_id", f"AF_{sp}"]], on="variant_id", how="outer")

    global_out = str(AF_OUT / f"freq_global{suffix}")
    if not Path(global_out + ".afreq").exists():
        run(
            [
                PLINK2,
                "--memory",
                str(plink_mem),
                "--pfile",
                str(pfile_prefix),
                "--extract",
                str(extract_file),
                "--freq",
                "--out",
                global_out,
            ],
            f"freq chr{chrom} global",
        )
    global_af = pd.read_csv(global_out + ".afreq", sep="\t")
    global_af = global_af.rename(columns={"#CHROM": "chr", "ID": "variant_id", "ALT_FREQS": "AF_global"})
    merged = merged.merge(global_af[["variant_id", "AF_global", "REF", "ALT"]], on="variant_id", how="left")

    superpop_path = AF_OUT / f"1000g_af_by_superpop{suffix}.parquet"
    merged.to_parquet(superpop_path, index=False)
    print(f"  Saved: {superpop_path.name} ({len(merged):,} variants)")
    maybe_write_legacy(chrom, superpop_path, "1000g_af_by_superpop.parquet")

    print("\nComputing AF divergence metrics ...")
    af_cols = [f"AF_{sp}" for sp in SUPERPOPS]
    af_mat = merged[af_cols].values.astype(float)
    pairs = list(combinations(range(len(SUPERPOPS)), 2))
    pairwise = np.array([np.abs(af_mat[:, i] - af_mat[:, j]) for i, j in pairs])
    merged["AF_max_diff"] = pairwise.max(axis=0)
    merged["AF_var"] = np.nanvar(af_mat, axis=1)

    p_bar = np.nanmean(af_mat, axis=1)
    var_p = np.nanvar(af_mat, axis=1)
    denom = p_bar * (1 - p_bar)
    fst = np.where(denom > 0, var_p / denom, 0.0)
    merged["FST_like"] = np.clip(fst, 0, 1)

    def ais(row: pd.Series) -> float:
        vals = np.array([row[c] for c in af_cols], dtype=float)
        vals = np.clip(vals, 1e-9, 1 - 1e-9)
        vals = vals / vals.sum()
        entropy = -np.sum(vals * np.log(vals))
        return float(entropy / np.log(len(vals)))

    merged["AIS"] = merged.apply(ais, axis=1)
    for (i, j), (sp1, sp2) in zip(pairs, [(SUPERPOPS[a], SUPERPOPS[b]) for a, b in pairs]):
        merged[f"AF_diff_{sp1}_{sp2}"] = np.abs(af_mat[:, i] - af_mat[:, j])

    divergence = merged.drop(columns=["OBS_CT"], errors="ignore")
    div_path = AF_OUT / f"af_divergence_features{suffix}.parquet"
    divergence.to_parquet(div_path, index=False)
    print(f"  Saved: {div_path.name} ({len(divergence):,} variants)")
    maybe_write_legacy(chrom, div_path, "af_divergence_features.parquet")

    n_labeled = len(all_variant_ids)
    n_with_af = divergence["AF_global"].notna().sum()
    pct = n_with_af / n_labeled * 100
    print(f"\nGate check chr{chrom}: {n_with_af:,}/{n_labeled:,} variants with AF features ({pct:.1f}%)")
    if pct < 80:
        print("  WARNING: <80% AF coverage. Investigate before proceeding.")
    else:
        print("  PASS: >80% AF coverage.")


if __name__ == "__main__":
    main()
