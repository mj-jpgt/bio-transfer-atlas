"""
FAIRGEN-Open Stage 4: Harmonize 9 PGS (GRCh38) → chr22 scoring pgen, run --score
===================================================================================
For each of the 9 selected PGS IDs:
  1. Read harmonized GRCh38 scoring file
  2. Filter to chr22
  3. Remove ambiguous A/T and C/G SNPs
  4. Match by chr:pos to chr22 scoring pvar
  5. Check effect allele concordance (direct / flip / mismatch)
  6. Write PLINK-ready harmonized score file
  7. Run plink2 --score
Outputs:
  data/processed/pgs_grch38/{PGS_ID}/{PGS_ID}.harmonized.tsv
  data/processed/pgs_grch38/{PGS_ID}/{PGS_ID}.match_report.tsv
  data/processed/pgs_grch38/harmonization_report.parquet
  data/processed/scores_grch38/{PGS_ID}.chr22.sscore
  data/processed/scores_grch38/score_matrix_grch38.parquet
"""
import subprocess
from pathlib import Path

import pandas as pd

root    = Path(__file__).resolve().parents[1]
plink2  = str(root / "tools/plink2/plink2.exe")
scores_raw = root / "data/raw/pgs_catalog/scores"
pgs_out    = root / "data/processed/pgs_grch38"
score_out  = root / "data/processed/scores_grch38"
pgs_out.mkdir(parents=True, exist_ok=True)
score_out.mkdir(parents=True, exist_ok=True)

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

AMBIGUOUS = {("A","T"),("T","A"),("C","G"),("G","C")}

# Load scoring pvar as chr:pos index
pvar_path = root / "data/interim/1000g_grch38/chr22.score.pvar.zst"
if not pvar_path.exists():
    pvar_path = root / "data/interim/1000g_grch38/chr22.score.pvar"

print("Loading GRCh38 chr22 scoring .pvar ...")
# PLINK2 pvar: header is '#CHROM  POS  ID  REF  ALT  ...'
pvar = pd.read_csv(pvar_path, sep="\t", comment="#", header=None,
                   usecols=[0,1,3,4],
                   names=["chr","pos","ref","alt"],
                   dtype={"chr": str, "pos": int, "ref": str, "alt": str})
pvar["chr"] = pvar["chr"].astype(str).str.replace("chr","", regex=False)
pvar_idx = pvar.set_index(["chr","pos"])
print(f"  scoring pvar variants: {len(pvar):,}")


def _find(cols, candidates):
    for c in candidates:
        if c in cols: return c
    return None


def harmonize_pgs(pgs_id: str, trait: str) -> dict:
    src = scores_raw / pgs_id / f"{pgs_id}_hmPOS_GRCh38.txt.gz"
    out_dir = pgs_out / pgs_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv    = out_dir / f"{pgs_id}.harmonized.tsv"
    report_tsv = out_dir / f"{pgs_id}.match_report.tsv"

    df = pd.read_csv(src, sep="\t", compression="gzip", comment="#",
                     dtype=str, low_memory=False)
    cols = df.columns.tolist()

    chr_col = _find(cols, ["hm_chr","chr_name"])
    pos_col = _find(cols, ["hm_pos","chr_position"])
    ea_col  = _find(cols, ["effect_allele"])
    oa_col  = _find(cols, ["other_allele"])
    ew_col  = _find(cols, ["effect_weight"])

    if not all([chr_col, pos_col, ea_col, ew_col]):
        missing = [n for n,v in [("chr",chr_col),("pos",pos_col),("ea",ea_col),("ew",ew_col)] if not v]
        return {"pgs_id": pgs_id, "error": f"missing columns: {missing}"}

    rename = {chr_col:"chr", pos_col:"pos", ea_col:"effect_allele", ew_col:"effect_weight"}
    if oa_col: rename[oa_col] = "other_allele"
    df = df[[c for c in rename]].rename(columns=rename)
    df["chr"] = df["chr"].astype(str).str.replace("chr","", regex=False)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["chr","pos","effect_allele","effect_weight"])
    df["pos"] = df["pos"].astype(int)

    n_original = len(df)
    df = df[df["chr"] == "22"].copy()
    n_chr22 = len(df)

    # Remove ambiguous SNPs
    n_ambig_removed = 0
    if "other_allele" in df.columns:
        pairs = list(zip(df["effect_allele"].str.upper(), df["other_allele"].fillna("").str.upper()))
        keep = [p not in AMBIGUOUS for p in pairs]
        n_ambig_removed = sum(not k for k in keep)
        df = df[keep].copy()
    n_unambig = len(df)

    # Match to pvar by chr:pos
    df = df.set_index(["chr","pos"])
    matched = df.join(pvar_idx[["ref","alt"]], how="inner").reset_index()
    n_pos_match = len(matched)

    # Allele concordance
    def allele_match(row):
        ea = row["effect_allele"].upper()
        ref = row["ref"].upper()
        alt = row["alt"].upper()
        if ea == alt: return "direct"
        if ea == ref: return "flip"
        return "mismatch"

    matched["_match"] = matched.apply(allele_match, axis=1)
    mismatches = (matched["_match"] == "mismatch").sum()
    matched = matched[matched["_match"] != "mismatch"].copy()
    n_allele_match = len(matched)

    # Build variant_id chr:pos:ref:alt and write PLINK score file
    matched["variant_id"] = (
        matched["chr"].astype(str) + ":" +
        matched["pos"].astype(str) + ":" +
        matched["ref"].astype(str) + ":" +
        matched["alt"].astype(str)
    )
    plink_df = matched[["variant_id","effect_allele","effect_weight"]].copy()
    plink_df.to_csv(out_tsv, sep="\t", index=False)

    # Match report
    report_df = matched[["variant_id","chr","pos","ref","alt","effect_allele","effect_weight","_match"]].copy()
    report_df.to_csv(report_tsv, sep="\t", index=False)

    rate_chr22 = n_allele_match / n_chr22 if n_chr22 > 0 else 0
    rate_global = n_allele_match / n_original if n_original > 0 else 0

    print(f"  {pgs_id} ({trait[:25]}) ... chr22={n_chr22} unambig={n_unambig} "
          f"matched={n_allele_match} rate={rate_chr22:.4f}")

    return {
        "pgs_id": pgs_id,
        "trait": trait,
        "n_variants_original": n_original,
        "n_chr22": n_chr22,
        "n_ambiguous_removed": n_ambig_removed,
        "n_pos_matched": n_pos_match,
        "n_allele_matched": n_allele_match,
        "n_mismatches_dropped": mismatches,
        "match_rate_chr22": round(rate_chr22, 4),
        "match_rate_global": round(rate_global, 6),
    }


def run_scoring(pgs_id: str) -> str | None:
    score_file = pgs_out / pgs_id / f"{pgs_id}.harmonized.tsv"
    if not score_file.exists():
        return None
    out_prefix = str(score_out / f"{pgs_id}.chr22")
    cmd = [
        plink2,
        "--pfile", str(root / "data/interim/1000g_grch38/chr22.score"),
        "--score", str(score_file), "1", "2", "3",
        "header", "cols=+scoresums",
        "--out", out_prefix,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  SCORE FAILED {pgs_id}: {r.stderr[-200:]}")
        return None
    sscore = Path(out_prefix + ".sscore")
    n = sum(1 for _ in open(sscore)) - 1 if sscore.exists() else 0
    print(f"  {pgs_id} scored → {n} samples")
    return str(sscore)


# ── main ────────────────────────────────────────────────────────────────────
print("\nHarmonizing 9 PGS (GRCh38) → chr22 ...")
reports = []
sscore_paths = {}
for pgs_id, trait in SELECTED_PGS.items():
    rep = harmonize_pgs(pgs_id, trait)
    reports.append(rep)

report_df = pd.DataFrame(reports)
report_df.to_parquet(pgs_out / "harmonization_report.parquet", index=False)
print("\nHarmonization report:")
print(report_df[["pgs_id","trait","n_chr22","n_allele_matched","match_rate_chr22"]].to_string(index=False))

# Gate: only score PGS with chr22 match rate >= 0.70
MIN_MATCH = 0.70
to_score = report_df[report_df["match_rate_chr22"] >= MIN_MATCH]
print(f"\nScoring {len(to_score)} PGS with match_rate_chr22 >= {MIN_MATCH} ...")
for _, row in to_score.iterrows():
    sscore = run_scoring(row["pgs_id"])
    if sscore:
        sscore_paths[row["pgs_id"]] = sscore

# Merge .sscore files
print("\nMerging .sscore → score_matrix_grch38.parquet ...")
panel = pd.read_parquet(root / "data/processed/sample_metadata_grch38.parquet")
merged = None
for pgs_id, path in sscore_paths.items():
    df = pd.read_csv(path, sep="\t", usecols=["#IID","SCORE1_SUM"])
    df = df.rename(columns={"#IID": "sample_id", "SCORE1_SUM": pgs_id})
    merged = df if merged is None else merged.merge(df, on="sample_id", how="outer")

merged = merged.merge(panel, on="sample_id", how="left")
out = score_out / "score_matrix_grch38.parquet"
merged.to_parquet(out, index=False)

score_cols = [c for c in merged.columns if c.startswith("PGS")]
print(f"\nscore_matrix_grch38.parquet: {len(merged)} samples × {len(score_cols)} scores")
print(merged[["sample_id","pop","super_pop"] + score_cols].head(5).to_string(index=False))
