#!/usr/bin/env bash
# Run remaining robustness artifact jobs on Lambda
set -euo pipefail
cd /lambda/nfs/geeg/fairness
source scripts/lambda_env.sh
mkdir -p results/logs/lambda results/tables
LOG=results/logs/lambda/robustness_artifacts.log
exec >>"$LOG" 2>&1
echo "=== START $(date -Is) ==="
PY=.venv/bin/python
$PY -u scripts/eval_trait_scale_portability.py || true
$PY -u scripts/run_external_validation_ci.py || true
$PY -u scripts/run_external_page_grch38.py --chroms 1,2,6,19,22 --size-per-chrom 15000 || true
$PY -u scripts/run_nested_ablation_and_paired_auroc.py || true
$PY -u scripts/eval_external_risk_ranking.py || true
$PY -u scripts/eval_sign_discordance_and_power.py || true
$PY -u scripts/eval_intervention_retention_controls.py || true
$PY -u scripts/run_popcorn_rg_available.py || true
echo "=== ROBUSTNESS_ARTIFACTS_DONE $(date -Is) ==="
