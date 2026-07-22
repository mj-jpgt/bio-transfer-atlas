"""
Phase 5: Run PLINK2 --score for each harmonized PGS on the scoring pgen (no MAF filter).
Outputs per score: data/processed/scores/<pgs_id>.chr22.sscore
Final merged output: data/processed/scores/score_matrix.parquet
"""
import subprocess
from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[1]
plink2 = str(root / "tools/plink2/plink2.exe")
pfile = str(root / "data/interim/1000g/chr22.score")
scores_dir = root / "data/processed/scores"
scores_dir.mkdir(parents=True, exist_ok=True)

pgs_dir = root / "data/processed/pgs"
harmonization_report = pd.read_parquet(pgs_dir / "harmonization_report.parquet")

# Only score files with a decent chr22 match rate
MIN_CHR22_MATCH = 0.70
to_score = harmonization_report[harmonization_report["match_rate_chr22"] >= MIN_CHR22_MATCH].copy()
print(f"Scores passing chr22 match rate >= {MIN_CHR22_MATCH}: {len(to_score)}")

results = []
for _, row in to_score.iterrows():
    pgs_id = row["pgs_id"]
    score_file = pgs_dir / pgs_id / f"{pgs_id}.harmonized.tsv"
    if not score_file.exists():
        print(f"  SKIP {pgs_id}: harmonized file missing")
        continue

    out_prefix = str(scores_dir / f"{pgs_id}.chr22")
    cmd = [
        plink2,
        "--pfile", pfile,
        "--score", str(score_file), "1", "2", "3",
        "header", "cols=+scoresums",
        "--out", out_prefix,
    ]
    print(f"  Scoring {pgs_id} ...", end=" ", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"FAILED\n{result.stderr[-300:]}")
        results.append({"pgs_id": pgs_id, "status": "failed"})
    else:
        sscore = Path(out_prefix + ".sscore")
        n_samples = sum(1 for _ in open(sscore)) - 1 if sscore.exists() else 0
        print(f"OK  ({n_samples} samples)")
        results.append({"pgs_id": pgs_id, "status": "ok", "sscore": str(sscore)})

# ── Merge all .sscore files into score_matrix.parquet ────────────────────────
print("\nMerging .sscore files into score_matrix.parquet ...")
scored = [r for r in results if r["status"] == "ok"]
if not scored:
    print("No scored files to merge.")
    exit(1)

# Load sample metadata for joining
panel = pd.read_parquet(root / "data/processed/sample_metadata.parquet")

merged = None
for r in scored:
    pgs_id = r["pgs_id"]
    sscore_path = r["sscore"]
    df = pd.read_csv(sscore_path, sep="\t", usecols=["#IID", "SCORE1_SUM"])
    df = df.rename(columns={"#IID": "sample_id", "SCORE1_SUM": pgs_id})
    if merged is None:
        merged = df
    else:
        merged = merged.merge(df, on="sample_id", how="outer")

merged = merged.merge(panel, on="sample_id", how="left")
out = root / "data/processed/scores/score_matrix.parquet"
merged.to_parquet(out, index=False)

score_cols = [c for c in merged.columns if c.startswith("PGS")]
print(f"score_matrix.parquet: {len(merged)} samples × {len(score_cols)} scores")
print(f"Score columns: {score_cols}")
print(merged[["sample_id", "pop", "super_pop"] + score_cols].head().to_string(index=False))
