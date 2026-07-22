#!/usr/bin/env bash
# After M4 apply finishes: re-run Duffy gate + reablate + figures + gate
set -euo pipefail
cd /lambda/nfs/geeg/fairness
source scripts/lambda_env.sh
PY=.venv/bin/python
exec > results/logs/lambda/m5_post_apply.log 2>&1
# wait until apply_intervention not running and score_intervention starts or done
for i in $(seq 1 120); do
  if ! pgrep -f 'apply_intervention.py' >/dev/null; then
    break
  fi
  sleep 30
done
$PY -u scripts/run_duffy_positive_control.py || true
# if score_intervention not already launched by m4 lane, launch
if ! pgrep -f 'score_intervention.py' >/dev/null && ! grep -q 'score int' results/logs/lambda/m4_lane_fast.log; then
  $PY -u scripts/score_intervention.py \
    --pgs-ids PGS000191,PGS004133,PGS001288 \
    --modes fst,maf,duffy_gate,random,filter_10 \
    --chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
    --jobs 6 --threads 2 --memory-mb 5000 --tag genomewide || true
fi
# wait for score_intervention
for i in $(seq 1 180); do
  if ! pgrep -f 'score_intervention.py' >/dev/null; then
    break
  fi
  sleep 60
done
$PY -u scripts/evaluate_intervention.py --tag genomewide || true
$PY -u scripts/summarize_mad_clades.py || true
$PY -u scripts/run_ldblock_and_baselines.py || true
$PY -u scripts/figures/make_paper_figures.py || true
$PY -u scripts/run_robustness_bootstrap.py || true
$PY -u scripts/gate_literature_roadmap.py || true
echo POST_APPLY_DONE
