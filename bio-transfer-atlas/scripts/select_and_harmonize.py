"""
Phase 3 + 4: Select best PGS per trait and harmonize to 1000G chr22 .pvar

NOTE: chr22 smoke test only. Not for genome-wide inference.
  - PCA variance explained values are inflated (single chromosome).
  - match_rate is chr22-only, not genome-wide.

Outputs:
  data/processed/pgs/score_selection_table.parquet
  data/processed/pgs/<pgs_id>/<pgs_id>.harmonized.tsv   (PLINK-ready)
  data/processed/pgs/<pgs_id>/<pgs_id>.match_report.tsv
  data/processed/pgs/harmonization_report.parquet
"""
from pathlib import Path

import pandas as pd

SMOKE_TEST_NOTE = "chr22 smoke test — not genome-wide"

root = Path(__file__).resolve().parents[1]
scores_raw = root / "data/raw/pgs_catalog/scores"
pgs_out = root / "data/processed/pgs"
pgs_out.mkdir(parents=True, exist_ok=True)

# ── 1. Load candidate table ───────────────────────────────────────────────────
cands = pd.read_csv(
    root / "data/raw/pgs_catalog/metadata/candidate_scores.tsv", sep="\t"
)
print(f"Total candidates: {len(cands)}")
print(cands["trait_query"].value_counts().to_string())

# ── 2. Score selection — explicit criteria (methods §5.2) ────────────────────
# Gate 1: harmonized file must exist on disk
def file_exists(pgs_id):
    p = scores_raw / pgs_id / f"{pgs_id}_hmPOS_GRCh37.txt.gz"
    return p.exists()

cands["file_exists"] = cands["pgs_id"].apply(file_exists)
cands = cands[cands["file_exists"]].copy()
print(f"\nWith harmonized file on disk: {len(cands)}")

# Gate 2: must have genome build info
cands = cands[cands["genome_build"].notna() & (cands["genome_build"] != "")].copy()
print(f"With genome_build metadata: {len(cands)}")

# Gate 3: must have n_variants parseable
cands["n_variants"] = pd.to_numeric(cands["n_variants_original"], errors="coerce")
cands = cands[cands["n_variants"].notna() & (cands["n_variants"] > 0)].copy()
print(f"With valid variant count: {len(cands)}")

# Scoring function — ranked by methods §5.2 priority:
# 1. ancestry diversity (multi/non-EUR > EUR-only)
# 2. metadata clarity (has development_ancestry)
# 3. variant count (more variants = broader genomic coverage, but penalise >10M imputation-heavy)
def selection_score(row):
    anc = str(row.get("development_ancestry", "")).lower()
    is_multi = int(any(k in anc for k in ["multi", "mae", "african", "east asian", "south asian", "admixed", "mao"]))
    is_non_eur_only = int("sas" in anc or "eas" in anc or "afr" in anc)
    has_anc_meta = int(anc not in ["", "nan", "none"])
    n = row["n_variants"]
    # penalise extremely large imputation-heavy scores slightly
    n_score = min(n, 5_000_000) / 5_000_000
    return (is_multi, is_non_eur_only, has_anc_meta, n_score)

cands["_sel_score"] = cands.apply(selection_score, axis=1)
cands = cands.sort_values("_sel_score", ascending=False)

# Take top 3 per trait — include: 1 multi-ancestry, 1 EUR, 1 non-EUR if available
selected = (
    cands.groupby("trait_query")
    .head(3)
    .reset_index(drop=True)
    .drop(columns=["_sel_score", "file_exists"])
)
selected["reason_selected"] = "ranked_by_ancestry_diversity_metadata_variant_count"
selected["scope"] = SMOKE_TEST_NOTE
selected.to_parquet(pgs_out / "score_selection_table.parquet", index=False)
print(f"\nSelected {len(selected)} scores across {selected['trait_query'].nunique()} traits:")
print(selected[["pgs_id", "trait_query", "n_variants", "development_ancestry"]].to_string(index=False))

# ── 3. Load SCORING pvar as reference (no MAF filter — preserves rare PGS variants) ──
print("\nLoading chr22 SCORING .pvar reference (no MAF filter)...")
pvar = pd.read_csv(
    root / "data/interim/1000g/chr22.score.pvar",
    sep="\t", comment="#",
    header=0,
    names=["chr", "pos", "id", "ref", "alt"],
    usecols=[0, 1, 2, 3, 4],
    dtype={"chr": str, "pos": int, "ref": str, "alt": str},
)
pvar["chr"] = pvar["chr"].str.replace("chr", "", regex=False)
pvar_idx = pvar.set_index(["chr", "pos"])
print(f"  scoring pvar variants: {len(pvar):,}  (includes MAF<0.01 variants)")

AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

def harmonize_pgs(pgs_id: str, trait: str) -> dict:
    src = scores_raw / pgs_id / f"{pgs_id}_hmPOS_GRCh37.txt.gz"
    out_dir = pgs_out / pgs_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / f"{pgs_id}.harmonized.tsv"
    report_tsv = out_dir / f"{pgs_id}.match_report.tsv"

    # Read PGS file, skip comment lines
    df = pd.read_csv(src, sep="\t", compression="gzip", comment="#",
                     dtype=str, low_memory=False)

    # Detect columns by name — prefer harmonized hm_ columns, fall back to raw
    cols = df.columns.tolist()
    def _find(candidates):
        for c in candidates:
            if c in cols:
                return c
        return None

    chr_col = _find(["hm_chr", "chr_name"])
    pos_col = _find(["hm_pos", "chr_position"])
    ea_col  = _find(["effect_allele"])
    oa_col  = _find(["other_allele"])
    ew_col  = _find(["effect_weight"])
    id_col  = _find(["hm_rsID", "rsID", "SNP"])

    if not all([chr_col, pos_col, ea_col, ew_col]):
        missing = [n for n, v in [("chr", chr_col), ("pos", pos_col), ("ea", ea_col), ("ew", ew_col)] if not v]
        return {"pgs_id": pgs_id, "error": f"missing required columns: {missing}"}

    # Build explicit rename map — no positional alignment
    rename_map = {chr_col: "chr", pos_col: "pos", ea_col: "effect_allele", ew_col: "effect_weight"}
    if oa_col:
        rename_map[oa_col] = "other_allele"
    keep_cols = [c for c in rename_map]
    df = df[keep_cols].rename(columns=rename_map)
    df["chr"] = df["chr"].astype(str).str.replace("chr", "", regex=False)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["chr", "pos", "effect_allele", "effect_weight"])
    df["pos"] = df["pos"].astype(int)

    n_original = len(df)

    # Filter to chr22 only
    df = df[df["chr"] == "22"].copy()
    n_chr22 = len(df)

    # Remove ambiguous SNPs
    if "other_allele" in df.columns:
        _pairs = list(zip(df["effect_allele"].str.upper(), df["other_allele"].fillna("").str.upper()))
        keep = [p not in AMBIGUOUS for p in _pairs]
        df = df[keep].copy()
    n_unambig = len(df)

    # Match to pvar by chr:pos
    df = df.set_index(["chr", "pos"])
    matched = df.join(pvar_idx[["ref", "alt"]], how="inner").reset_index()
    n_pos_match = len(matched)

    # Check allele concordance; flip if needed
    def allele_match(row):
        ea = row["effect_allele"].upper()
        oa = row.get("other_allele", "")
        ref = row["ref"].upper() if pd.notna(row["ref"]) else ""
        alt = row["alt"].upper() if pd.notna(row["alt"]) else ""
        if ea == alt:
            return "direct"
        if ea == ref:
            return "flip"  # effect allele is REF — unusual but keep
        return "mismatch"

    matched["_match"] = matched.apply(allele_match, axis=1)
    matched = matched[matched["_match"] != "mismatch"].copy()
    n_allele_match = len(matched)

    # Build PLINK score file: variant_id effect_allele effect_weight
    # variant_id must match the IDs we set: chr:pos:ref:alt
    matched["variant_id"] = (
        matched["chr"].astype(str) + ":" +
        matched["pos"].astype(str) + ":" +
        matched["ref"] + ":" + matched["alt"]
    )
    plink_df = matched[["variant_id", "effect_allele", "effect_weight"]].copy()
    plink_df.to_csv(out_tsv, sep="\t", index=False)

    report = {
        "pgs_id": pgs_id,
        "trait": trait,
        "n_variants_original": n_original,
        "n_variants_chr22": n_chr22,
        "n_unambiguous": n_unambig,
        "n_position_matched": n_pos_match,
        "n_allele_matched": n_allele_match,
        "match_rate": round(n_allele_match / n_original, 4) if n_original else 0,
        "match_rate_chr22": round(n_allele_match / n_chr22, 4) if n_chr22 else 0,
    }
    pd.DataFrame([report]).to_csv(report_tsv, sep="\t", index=False)
    return report

# ── 4. Harmonize each selected score ─────────────────────────────────────────
print("\nHarmonizing scores...")
reports = []
for _, row in selected.iterrows():
    print(f"  {row['pgs_id']} ({row['trait_query']}) ...", end=" ")
    rep = harmonize_pgs(row["pgs_id"], row["trait_query"])
    print(f"chr22_matched={rep.get('n_allele_matched','ERR')}  rate={rep.get('match_rate_chr22','ERR')}")
    reports.append(rep)

match_df = pd.DataFrame(reports)
match_df.to_parquet(pgs_out / "harmonization_report.parquet", index=False)
print("\nHarmonization report:")
print(match_df[["pgs_id", "trait", "n_variants_original", "n_allele_matched", "match_rate", "match_rate_chr22"]].to_string(index=False))
