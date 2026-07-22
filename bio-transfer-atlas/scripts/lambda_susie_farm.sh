#!/usr/bin/env bash
# Farm SuSiE across traits x ancestries for one chrom.
set -euo pipefail
ROOT=/lambda/nfs/geeg/fairness
cd "$ROOT"
source .venv/bin/activate
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
CHROM="${1:-22}"
LOG=results/logs/lambda/m3b_susie_farm.log
{
  for trait in T2D CAD BMI LDL; do
    for anc in EUR AFR; do
      echo "=== ${trait} ${anc} chr${CHROM} $(date -Is) ==="
      python scripts/run_polyfun_susie.py --trait "$trait" --chrom "$CHROM" --anc "$anc" --jobs 12 --max-blocks 100
    done
  done
  echo "=== build tiers ==="
  python scripts/build_finemap_tier_labels.py --tag genomewide
  echo "=== eval tiers ==="
  python scripts/eval_finemap_tiers_lean.py || python scripts/eval_finemap_tiers.py || true
  echo SUSIE_FARM_DONE
} 2>&1 | tee -a "$LOG"
