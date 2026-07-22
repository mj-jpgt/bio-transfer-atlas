#!/usr/bin/env bash
set -euo pipefail
ROOT=/lambda/nfs/geeg/fairness
cd "$ROOT"
# shellcheck disable=SC1091
source scripts/lambda_env.sh
mkdir -p results/logs/lambda

pkill -f 'scripts/evaluate_intervention.py' 2>/dev/null || true
sleep 2

nohup .venv/bin/python -u scripts/evaluate_intervention.py \
  --baseline data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet \
  --labels data/modeling/master_variant_table_genomewide_genomewide.parquet \
  --intervention-root data/processed/pgs_grch38_intervention_genomewide \
  --intervention-scores data/processed/scores_grch38_intervention_genomewide \
  --tag genomewide \
  --n-boot 100 \
  --out-csv results/tables/intervention_results.genomewide.csv \
  --out-summary results/tables/intervention_summary.genomewide.txt \
  > results/logs/lambda/m2_evaluate_intervention.log 2>&1 &
echo "EVAL_PID=$!"
sleep 20
tail -40 results/logs/lambda/m2_evaluate_intervention.log || true
