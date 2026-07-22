#!/usr/bin/env bash
# Restart SuSiE with fixed --r2-unphased real LD (priority chroms)
set -euo pipefail
cd /lambda/nfs/geeg/fairness
source scripts/lambda_env.sh
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
# stop old farm
pkill -f 'run_polyfun_susie.py' || true
pkill -f 'lambda_m5_susie.sh' || true
sleep 2
exec > results/logs/lambda/m5_susie_reald.log 2>&1
PY=.venv/bin/python
for trait in T2D CAD BMI LDL; do
  for anc in EUR AFR; do
    for chrom in 22 1 6 19; do
      echo "=== REAL_LD $trait $anc chr$chrom ==="
      $PY -u scripts/run_polyfun_susie.py \
        --trait "$trait" --chrom "$chrom" --anc "$anc" \
        --jobs 4 --max-blocks 25 --prefer-real-ld || true
    done
  done
done
$PY -u scripts/_postprocess_susie_cs.py || true
$PY -u scripts/build_finemap_tier_labels.py --tag genomewide_susie || true
$PY -u scripts/eval_finemap_tiers_lean.py || true
if test -f data/labels/susie/susie_real_ld_status_summary.csv; then
  cp data/labels/susie/susie_real_ld_status_summary.csv results/tables/
fi
echo SUSIE_REAL_LD_DONE
