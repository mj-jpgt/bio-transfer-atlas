#!/usr/bin/env bash
set -euo pipefail
ROOT=/lambda/nfs/geeg/fairness
cd "$ROOT"
# shellcheck disable=SC1091
source scripts/lambda_env.sh
mkdir -p results/logs/lambda

# Stop any leftover sequential scorers / orphan plink
pkill -f 'scripts/score_intervention.py' 2>/dev/null || true
pkill -f 'scripts/run_genomewide_downstream.py --step atlas,intervention' 2>/dev/null || true
# orphan plink from old score run
pkill -f 'scores_grch38_intervention_genomewide/_sscore_work' 2>/dev/null || true
sleep 2

N_DONE=$(ls data/processed/scores_grch38_intervention_genomewide/_sscore_work/*.sscore 2>/dev/null | wc -l || echo 0)
echo "[launch] existing sscores=${N_DONE}"

nohup .venv/bin/python -u scripts/score_intervention.py \
  --chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
  --intervention-root data/processed/pgs_grch38_intervention_genomewide \
  --out-dir data/processed/scores_grch38_intervention_genomewide \
  --memory-mb 6000 --jobs 12 --threads 3 --tag genomewide \
  > results/logs/lambda/m2_score_intervention.log 2>&1 &
echo "SCORE_PID=$!"
sleep 15
pgrep -af 'score_intervention|plink2' | grep -v grep | head -20 || true
tail -30 results/logs/lambda/m2_score_intervention.log || true
free -h | head -2
