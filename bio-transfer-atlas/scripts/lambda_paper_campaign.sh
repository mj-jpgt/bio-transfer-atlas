#!/usr/bin/env bash
# Paper completion campaign on Lambda: M4 scoring + M5 subpop/rg/susie/gnn lanes.
set -euo pipefail
ROOT=/lambda/nfs/geeg/fairness
cd "$ROOT"
# shellcheck disable=SC1091
source scripts/lambda_env.sh
mkdir -p results/logs/lambda

echo "[campaign] sync check"
.venv/bin/python - <<'PY'
from pathlib import Path
for p in ["PGS000191","PGS004133","PGS001288"]:
    f=Path(f"data/raw/pgs_catalog/scores/{p}/{p}_hmPOS_GRCh38.txt.gz")
    print(p, "raw", f.exists(), f.stat().st_size if f.exists() else 0)
    h=Path(f"data/processed/pgs_grch38/{p}/{p}.harmonized.tsv")
    print(p, "harm", h.exists())
PY

# --- Lane A: harmonize + score expansion PGS (CPU) ---
if ! test -f data/processed/pgs_grch38/PGS000191/PGS000191.harmonized.tsv; then
  echo "[campaign] harmonize expansion PGS"
  nohup .venv/bin/python -u scripts/harmonize_pgs_genomewide.py \
    > results/logs/lambda/m4_harmonize.log 2>&1 &
  echo HARMONIZE_PID=$!
  wait $!
fi

echo "[campaign] score expansion PGS into genomewide matrix (resume-friendly)"
nohup .venv/bin/python -u scripts/score_genomewide.py \
  --pgs-ids PGS000191,PGS004133,PGS001288 \
  --chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
  --jobs 8 --threads 3 --memory-mb 6000 \
  --out data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet \
  > results/logs/lambda/m4_score_expand.log 2>&1 &
echo SCORE_PID=$!

# --- Lane B: subpop AF chr1-7 ---
nohup .venv/bin/python -u scripts/compute_subpop_af_features.py \
  --chroms 1,2,3,4,5,6,7 --jobs 4 --memory-mb 8192 \
  > results/logs/lambda/m5_subpop_chr1_7.log 2>&1 &
echo SUBPOP_PID=$!

# --- Lane C: Popcorn autosomal (long) ---
nohup .venv/bin/python -u scripts/run_popcorn_rg.py --autosome \
  > results/logs/lambda/m5_popcorn_autosome.log 2>&1 &
echo POPCORN_PID=$!

echo "[campaign] launched SCORE=$SCORE_PID SUBPOP=$SUBPOP_PID POPCORN=$POPCORN_PID"
sleep 5
pgrep -af 'score_genomewide|compute_subpop|run_popcorn' | grep -v grep | head -20 || true
