#!/usr/bin/env bash
# SuSiE farm with PRIMARY signed LD (--r-unphased).
set -euo pipefail
cd /lambda/nfs/geeg/fairness
source scripts/lambda_env.sh
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
mkdir -p results/logs/lambda data/labels/susie
LOG=results/logs/lambda/m5_susie_signed_ld.log
exec >>"$LOG" 2>&1
echo "=== START $(date -Is) ==="
pkill -f 'run_polyfun_susie.py' || true
sleep 1
PY=.venv/bin/python
$PY - <<'PY'
import json
from pathlib import Path
p = Path("data/labels/susie/susie_primary_meta.json")
p.parent.mkdir(parents=True, exist_ok=True)
p.write_text(json.dumps({
    "primary_requires": "signed_ld",
    "legacy_r2_unsigned_abs_r": "non_primary",
    "identity_ld": "non_primary",
    "absz_weight": "non_primary",
}, indent=2), encoding="utf-8")
print("wrote", p)
PY
for trait in T2D CAD BMI LDL; do
  for anc in EUR AFR; do
    for chrom in 22 1 6 19; do
      echo "=== SIGNED_LD $trait $anc chr$chrom $(date -Is) ==="
      $PY -u scripts/run_polyfun_susie.py \
        --trait "$trait" --chrom "$chrom" --anc "$anc" \
        --jobs 4 --max-blocks 25 --prefer-real-ld || true
    done
  done
done
$PY -u scripts/_postprocess_susie_cs.py || true
$PY -u scripts/build_finemap_tier_labels.py --tag genomewide_susie || true
$PY -u scripts/build_finemap_tier_labels_zlead.py --tag genomewide_susie_zlead || true
$PY -u scripts/eval_finemap_tiers_lean.py || true
if test -f data/labels/susie/susie_real_ld_status_summary.csv; then
  cp data/labels/susie/susie_real_ld_status_summary.csv results/tables/
fi
echo "=== SUSIE_SIGNED_LD_DONE $(date -Is) ==="
