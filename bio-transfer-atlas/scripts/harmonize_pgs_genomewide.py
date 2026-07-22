"""
Harmonize 9 PGS scores against the full genome-wide 1000G GRCh38 scoring pvars.

Memory-efficient strategy:
  1. Load all 9 PGS files upfront (small: ~1 M rows each).
  2. For each chromosome, do ONE chunked scan of chrN.score.pvar, keeping
     only rows whose position appears in any PGS file for that chromosome.
  3. Match each PGS file against the filtered pvar subset.

Peak RSS stays well under 4 GB; total pvar I/O is O(22 chromosomes) not
O(9 PGS files × 22 chromosomes).

Produces: data/processed/pgs_grch38/{PGS_ID}/{PGS_ID}.harmonized.tsv
Run once before run_genomewide.py.  Safe to re-run.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCORES_RAW = ROOT / "data/raw/pgs_catalog/scores"
PGS_OUT = ROOT / "data/processed/pgs_grch38"
INTERIM = ROOT / "data/interim/1000g_grch38"
PGS_OUT.mkdir(parents=True, exist_ok=True)

SELECTED_PGS = {
    "PGS000018": "type 2 diabetes",
    "PGS004696": "coronary artery disease",
    "PGS004698": "coronary artery disease",
    "PGS003897": "body mass index",
    "PGS002853": "LDL cholesterol",
    "PGS002858": "LDL cholesterol",
    "PGS003092": "body mass index",
    "PGS000014": "coronary artery disease",
    "PGS004840": "type 2 diabetes",
    "PGS000191": "white blood cell count",
    "PGS004133": "rheumatoid arthritis",
    "PGS001288": "inflammatory bowel disease",
}

AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def _find(cols: list[str], candidates: list[str]) -> str | None:
    for c in candidates:
        if c in cols:
            return c
    return None


def available_chroms() -> list[str]:
    chroms = [
        p.stem.split(".")[0].replace("chr", "")
        for p in sorted(INTERIM.glob("chr*.score.pvar"))
    ]
    return sorted((c for c in chroms if c.isdigit()), key=int)


def read_and_clean_pgs(pgs_id: str) -> pd.DataFrame | None:
    """Load one PGS file; return cleaned DataFrame or None on error."""
    src = SCORES_RAW / pgs_id / f"{pgs_id}_hmPOS_GRCh38.txt.gz"
    if not src.exists():
        print(f"  {pgs_id}: source missing, skipping")
        return None
    df = pd.read_csv(src, sep="\t", compression="gzip", comment="#", dtype=str, low_memory=False)
    cols = df.columns.tolist()
    chr_col = _find(cols, ["hm_chr", "chr_name"])
    pos_col = _find(cols, ["hm_pos", "chr_position"])
    ea_col  = _find(cols, ["effect_allele"])
    oa_col  = _find(cols, ["other_allele"])
    ew_col  = _find(cols, ["effect_weight"])
    if not all([chr_col, pos_col, ea_col, ew_col]):
        print(f"  {pgs_id}: missing required columns, skipping")
        return None
    rename = {chr_col: "chr", pos_col: "pos", ea_col: "effect_allele", ew_col: "effect_weight"}
    if oa_col:
        rename[oa_col] = "other_allele"
    df = df[[c for c in rename]].rename(columns=rename)
    df["chr"] = df["chr"].astype(str).str.replace("chr", "", regex=False)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["chr", "pos", "effect_allele", "effect_weight"]).copy()
    df["pos"] = df["pos"].astype(int)
    # Keep autosomes only; remove ambiguous SNPs
    df = df[df["chr"].str.match(r"^\d+$")].copy()
    if "other_allele" in df.columns:
        pairs = list(zip(df["effect_allele"].str.upper(), df["other_allele"].fillna("").str.upper()))
        df = df[[p not in AMBIGUOUS for p in pairs]].copy()
    return df


def pvar_cache_path(chrom: str) -> Path:
    return INTERIM / f"chr{chrom}.score_lookup.parquet"


def build_pvar_cache(chrom: str) -> None:
    """Convert chrN.score.pvar → compact parquet (chr/pos/ref/alt only).

    Reads the large pvar text file once and caches it as a ~20 MB parquet,
    making all subsequent lookups instant and memory-safe.
    """
    cache = pvar_cache_path(chrom)
    if cache.exists():
        return
    pvar_path = INTERIM / f"chr{chrom}.score.pvar"
    print(f"  chr{chrom}: building pvar cache (reading {pvar_path.stat().st_size / 1e9:.1f} GB pvar) ...", flush=True)
    parts: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        pvar_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 3, 4],
        names=["chr", "pos", "ref", "alt"],
        dtype={"chr": str, "pos": int, "ref": str, "alt": str},
        chunksize=1_000_000,
    ):
        chunk["chr"] = chunk["chr"].astype(str).str.replace("chr", "", regex=False)
        parts.append(chunk)
    pvar = pd.concat(parts, ignore_index=True)
    pvar.to_parquet(cache, index=False)
    print(f"    cached {len(pvar):,} variants -> {cache.name} ({cache.stat().st_size / 1e6:.0f} MB)", flush=True)


def scan_pvar_for_positions(chrom: str, pos_set: set[int]) -> pd.DataFrame:
    """Load compact pvar cache and filter to positions in pos_set."""
    cache = pvar_cache_path(chrom)
    if not cache.exists():
        build_pvar_cache(chrom)
    pvar = pd.read_parquet(cache)
    return pvar[pvar["pos"].isin(pos_set)].copy()


def match_pgs_to_pvar(pgs_sub: pd.DataFrame, pvar_sub: pd.DataFrame) -> pd.DataFrame:
    pvar_idx = pvar_sub.set_index(["chr", "pos"])
    matched = pgs_sub.set_index(["chr", "pos"]).join(pvar_idx[["ref", "alt"]], how="inner").reset_index()

    def _match(row: pd.Series) -> str:
        ea = row["effect_allele"].upper()
        if ea == str(row["alt"]).upper():
            return "direct"
        if ea == str(row["ref"]).upper():
            return "flip"
        return "mismatch"

    matched["_m"] = matched.apply(_match, axis=1)
    matched = matched[matched["_m"] != "mismatch"].copy()
    matched["variant_id"] = (
        matched["chr"].astype(str) + ":"
        + matched["pos"].astype(str) + ":"
        + matched["ref"].astype(str) + ":"
        + matched["alt"].astype(str)
    )
    return matched[["variant_id", "effect_allele", "effect_weight"]].copy()


def already_done() -> bool:
    """Return True if all harmonized TSVs already exist (skip re-run)."""
    report = PGS_OUT / "harmonization_report_genomewide.parquet"
    if not report.exists():
        return False
    return all((PGS_OUT / pid / f"{pid}.harmonized.tsv").exists() for pid in SELECTED_PGS)


def main() -> None:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--force", action="store_true", help="Re-run even if TSVs already exist")
    args = p.parse_args()

    if not args.force and already_done():
        print("Genome-wide harmonized TSVs already present — skipping (use --force to re-run).")
        return

    chroms = available_chroms()
    if not chroms:
        raise FileNotFoundError(
            f"No chr*.score.pvar files found in {INTERIM}. "
            "Run preprocess_grch38_genomewide.py first."
        )
    print(f"Available chromosomes: {chroms}\n")

    # --- Step 1: build pvar caches for all chromosomes (one-time, reused below) ---
    print("Building pvar caches ...")
    for chrom in chroms:
        if not (INTERIM / f"chr{chrom}.score.pvar").exists():
            print(f"  chr{chrom}: no pvar, skipping")
            continue
        build_pvar_cache(chrom)

    # --- Step 2: process one PGS file at a time to cap memory usage ---
    # With caches built, each pvar load is a small parquet read (~20 MB).
    matched_by_pgs: dict[str, list[pd.DataFrame]] = {pid: [] for pid in SELECTED_PGS}
    n_original: dict[str, int] = {}

    for pgs_id, trait in SELECTED_PGS.items():
        print(f"\nProcessing {pgs_id} ({trait}) ...", flush=True)
        pgs_df = read_and_clean_pgs(pgs_id)
        if pgs_df is None:
            n_original[pgs_id] = 0
            continue
        n_original[pgs_id] = len(pgs_df)
        print(f"  {len(pgs_df):,} variants loaded", flush=True)

        for chrom in chroms:
            pgs_sub = pgs_df[pgs_df["chr"] == chrom]
            if pgs_sub.empty:
                continue
            cache = pvar_cache_path(chrom)
            if not cache.exists():
                continue
            pvar = pd.read_parquet(cache)
            pos_set = set(pgs_sub["pos"].tolist())
            pvar_sub = pvar[pvar["pos"].isin(pos_set)]
            del pvar  # free immediately
            if pvar_sub.empty:
                continue
            hits = match_pgs_to_pvar(pgs_sub, pvar_sub)
            if not hits.empty:
                matched_by_pgs[pgs_id].append(hits)

        del pgs_df  # free PGS data before loading next one

    # --- Step 3: write harmonized TSVs ---
    print("\nWriting harmonized TSVs ...")
    reports = []
    for pgs_id, trait in SELECTED_PGS.items():
        out_dir = PGS_OUT / pgs_id
        out_dir.mkdir(parents=True, exist_ok=True)
        out_tsv = out_dir / f"{pgs_id}.harmonized.tsv"
        frames = matched_by_pgs.get(pgs_id, [])
        result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
            columns=["variant_id", "effect_allele", "effect_weight"]
        )
        result.to_csv(out_tsv, sep="\t", index=False)
        n_matched = len(result)
        n_orig = n_original.get(pgs_id, 0)
        rate = n_matched / n_orig if n_orig else 0
        print(f"  {pgs_id} ({trait[:30]}) ... {n_matched:,}/{n_orig:,}  rate={rate:.4f}")
        reports.append({
            "pgs_id": pgs_id, "trait": trait,
            "n_variants_original": n_orig,
            "n_allele_matched": n_matched,
            "match_rate_global": round(rate, 6),
        })

    pd.DataFrame(reports).to_parquet(
        PGS_OUT / "harmonization_report_genomewide.parquet", index=False
    )
    print("\nDone.")


if __name__ == "__main__":
    main()
