#!/usr/bin/env bash
# Re-farm remaining traits with z-fallback SuSiE, then eval.
set -euo pipefail
cd /lambda/nfs/geeg/fairness
source .venv/bin/activate
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
for trait in T2D CAD BMI LDL; do
  for anc in EUR AFR; do
    echo "=== $trait $anc ==="
    python scripts/run_polyfun_susie.py --trait "$trait" --chrom 22 --anc "$anc" --jobs 16 --max-blocks 80
  done
done
python scripts/_postprocess_susie_cs.py
python scripts/eval_finemap_tiers_lean.py
echo FARM2_DONE
