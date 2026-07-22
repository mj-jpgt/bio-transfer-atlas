#!/usr/bin/env bash
# Master Lambda campaign: overlapping waves for M1-M3 + GPU lane.
set -euo pipefail
ROOT="${BTA_ROOT:-/lambda/nfs/geeg/fairness}"
cd "$ROOT"
chmod +x scripts/lambda_*.sh 2>/dev/null || true
# shellcheck disable=SC1091
source scripts/lambda_env.sh
LOG="$ROOT/results/logs/lambda"
mkdir -p "$LOG" data/annotations data/raw/1000g/vcf_grch38 data/interim/1000g_grch38 data/features/baselines
PY="${PYTHON:-python}"

echo "=== BOOTSTRAP (if needed) ==="
if [[ ! -x tools/plink2/plink2 ]] || ! command -v Rscript >/dev/null 2>&1; then
  bash scripts/lambda_bootstrap.sh 2>&1 | tee "$LOG/bootstrap.log"
fi
# shellcheck disable=SC1091
source scripts/lambda_env.sh

echo "=== WAVE overlap: M1 AM + M2 VCF download ==="
nohup $PY scripts/join_alphamissense.py \
  --variant-list data/modeling/_tmp_ldblock_associated_sample.parquet \
  --out data/annotations/alphamissense_grch38.parquet \
  --min-free-gb 0 >"$LOG/m1_am.log" 2>&1 &
echo $! >"$LOG/m1_am.pid"

nohup $PY scripts/download_1000g_grch38_all_chrs.py --chroms 1-7 --jobs 4 \
  >"$LOG/m2_vcf.log" 2>&1 &
echo $! >"$LOG/m2_vcf.pid"

echo "Waiting for AlphaMissense join ..."
wait "$(cat "$LOG/m1_am.pid")" || true
$PY scripts/eval_vep_af_interaction.py --min-free-gb 0 2>/dev/null \
  || $PY scripts/eval_vep_af_interaction.py >"$LOG/m1_eval.log" 2>&1

echo "Waiting for VCF downloads ..."
wait "$(cat "$LOG/m2_vcf.pid")" || true

echo "=== M2 rebuild (2 parallel) ==="
$PY scripts/rebuild_chr1_7_score_pgens.py \
  --chroms 1,2,3,4,5,6,7 \
  --memory-mb "${BTA_PLINK_MEMORY_MB}" \
  --threads "${BTA_PLINK_THREADS}" \
  --jobs "${BTA_REBUILD_JOBS}" \
  2>&1 | tee "$LOG/m2_rebuild.log"

echo "=== GPU SHAP (background) while M3 starts ==="
nohup $PY scripts/run_shap_attribution.py >"$LOG/gpu_shap.log" 2>&1 &
echo $! >"$LOG/gpu_shap.pid"

echo "=== M3A Popcorn/rg ==="
$PY scripts/run_popcorn_rg.py --chroms 22 --traits T2D,CAD,BMI,LDL 2>&1 | tee "$LOG/m3a_popcorn.log"

echo "=== M3B SuSiE farm (chr22, 4 traits) ==="
for trait in T2D CAD BMI LDL; do
  for anc in EUR AFR; do
    $PY scripts/run_polyfun_susie.py --trait "$trait" --chrom 22 --anc "$anc" --jobs 16 --max-blocks 150 \
      2>&1 | tee -a "$LOG/m3b_susie.log"
  done
done
$PY scripts/build_finemap_tier_labels.py --tag genomewide 2>&1 | tee -a "$LOG/m3b_susie.log"
$PY scripts/eval_finemap_tiers_lean.py 2>&1 | tee -a "$LOG/m3b_tiers_eval.log" || \
  $PY scripts/eval_finemap_tiers.py 2>&1 | tee -a "$LOG/m3b_tiers_eval.log" || true

echo "=== M3A ablation with RG_REAL ==="
$PY scripts/run_ldblock_and_baselines.py 2>&1 | tee "$LOG/m3a_ablation.log"

echo "=== M3C subpop AF ==="
$PY scripts/compute_subpop_af_features.py --chroms 22 --memory-mb 8192 2>&1 | tee "$LOG/m3c_subpop.log" || true
# scale if score pgens exist for more chroms
HAVE=$(ls data/interim/1000g_grch38/chr*.score.pgen 2>/dev/null | wc -l)
if [[ "$HAVE" -ge 8 ]]; then
  $PY scripts/compute_subpop_af_features.py --chroms 1-22 --jobs 4 --memory-mb 8192 \
    2>&1 | tee -a "$LOG/m3c_subpop.log" || true
fi

echo "=== M2 downstream score/atlas/intervention if chr1-7 ok ==="
if [[ -f results/tables/score_pgen_chr1_7_rebuild_status.csv ]]; then
  if grep -qv 'blocked\|failed' results/tables/score_pgen_chr1_7_rebuild_status.csv; then
    $PY scripts/run_genomewide_downstream.py --step score,atlas,intervention \
      --tag genomewide --chroms 1-22 --score-chroms 1-22 \
      2>&1 | tee "$LOG/m2_downstream.log" || true
  fi
fi

wait "$(cat "$LOG/gpu_shap.pid")" 2>/dev/null || true
echo "CAMPAIGN_DONE"
date
