#!/usr/bin/env bash
# Remaining robustness Lambda jobs: Duffy ACKR1, matched MC controls, PAGE QC
set -euo pipefail
ROOT=/lambda/nfs/geeg/fairness
cd "$ROOT"
export PYTHONPATH="$ROOT/scripts:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1
LOG="$ROOT/results/logs/lambda/remaining_robustness_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")" "$ROOT/results/tables"
exec > >(tee -a "$LOG") 2>&1
echo "START remaining robustness $(date -Is)"
echo "LOG=$LOG"

PY=python3
if [[ -x "$ROOT/.venv/bin/python" ]]; then PY="$ROOT/.venv/bin/python"; fi
if [[ -x "$ROOT/venv/bin/python" ]]; then PY="$ROOT/venv/bin/python"; fi

echo "=== 1) Duffy homozygous + ACKR1 decomposition ==="
$PY scripts/run_duffy_positive_control.py || echo "DUFFY_FAILED=$?"

echo "=== 2) Intervention matched-n/mass Monte Carlo + LOSO ==="
$PY scripts/run_intervention_matched_random_controls.py --n-draws 500 || echo "MATCHED_FAILED=$?"

echo "=== 3) PAGE allele QC ladder ==="
$PY scripts/run_external_page_grch38.py || echo "PAGE_FAILED=$?"

echo "=== ARTIFACT CHECK ==="
ls -la results/tables/duffy_ackr1_score_decomposition.csv \
  results/tables/duffy_positive_control_genomewide.csv \
  results/tables/intervention_matched_random_controls.csv \
  results/tables/intervention_loso_mad_by_mode.csv \
  results/tables/external_page_validation.csv \
  results/tables/external_page_qc_counts.csv 2>&1 || true

echo "DONE remaining robustness $(date -Is)"
