#!/usr/bin/env bash
set -euo pipefail
cd /lambda/nfs/geeg/fairness
source scripts/lambda_env.sh
PY=.venv/bin/python
mkdir -p results/logs/lambda
exec > results/logs/lambda/m4_lane_fast.log 2>&1
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] lean harmonize"
$PY -u scripts/harmonize_pgs_lean.py
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] score expand"
$PY -u scripts/score_genomewide.py \
  --pgs-ids PGS000191,PGS004133,PGS001288 \
  --chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
  --jobs 8 --threads 3 --memory-mb 6000 \
  --out data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] duffy"
$PY -u scripts/run_duffy_positive_control.py || true
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] apply int"
$PY -u scripts/apply_intervention.py --modes filter_10,random,fst,maf,duffy_gate || true
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] score int"
$PY -u scripts/score_intervention.py \
  --pgs-ids PGS000191,PGS004133,PGS001288 \
  --modes fst,maf,duffy_gate,random,filter_10 \
  --chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
  --jobs 6 --threads 2 --memory-mb 5000 --tag genomewide || true
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] eval"
$PY -u scripts/evaluate_intervention.py --tag genomewide || true
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] mad clades"
$PY -u scripts/summarize_mad_clades.py || true
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] M4_LANE_DONE"
