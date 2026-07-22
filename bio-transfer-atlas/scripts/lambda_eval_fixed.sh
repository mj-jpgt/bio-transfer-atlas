#!/usr/bin/env bash
cd /lambda/nfs/geeg/fairness
source scripts/lambda_env.sh
.venv/bin/python -u scripts/evaluate_intervention.py --tag genomewide \
  --baseline data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet \
  > results/logs/lambda/m4_eval_fixed.log 2>&1
.venv/bin/python -u scripts/summarize_mad_clades.py >> results/logs/lambda/m4_eval_fixed.log 2>&1
.venv/bin/python -u scripts/run_duffy_positive_control.py >> results/logs/lambda/m4_eval_fixed.log 2>&1
echo EVAL_FIXED_DONE >> results/logs/lambda/m4_eval_fixed.log
