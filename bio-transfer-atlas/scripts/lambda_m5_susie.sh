#!/usr/bin/env bash
# Susie real-LD priority chroms + postprocess
set -euo pipefail
cd /lambda/nfs/geeg/fairness
source scripts/lambda_env.sh
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
exec > results/logs/lambda/m5_susie_farm.log 2>&1
PY=.venv/bin/python
for trait in T2D CAD BMI LDL; do
  for anc in EUR AFR; do
    for chrom in 1 2 3 6 7 19 22; do
      echo "=== $trait $anc chr$chrom ==="
      $PY -u scripts/run_polyfun_susie.py \
        --trait "$trait" --chrom "$chrom" --anc "$anc" \
        --jobs 6 --max-blocks 30 --prefer-real-ld || true
    done
  done
done
$PY -u scripts/_postprocess_susie_cs.py || true
$PY -u scripts/build_finemap_tier_labels.py --tag genomewide_susie || true
$PY -u scripts/eval_finemap_tiers_lean.py || true
# copy status summary to results/tables
if test -f data/labels/susie/susie_real_ld_status_summary.csv; then
  cp data/labels/susie/susie_real_ld_status_summary.csv results/tables/
fi
echo SUSIE_FARM_DONE
