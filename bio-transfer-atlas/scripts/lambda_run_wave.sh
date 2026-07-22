#!/usr/bin/env bash
# Launch overlapping Lambda waves. Logs -> results/logs/lambda/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck disable=SC1091
source "$ROOT/scripts/lambda_env.sh"
WAVE="${1:-help}"
LOG="$BTA_ROOT/results/logs/lambda"
PY="${PYTHON:-python3}"

run_bg() {
  local name="$1"; shift
  echo "[launch] $name -> $LOG/${name}.log"
  nohup "$@" >"$LOG/${name}.log" 2>&1 &
  echo $! >"$LOG/${name}.pid"
  echo "  pid=$(cat "$LOG/${name}.pid")"
}

case "$WAVE" in
  wave0)
    echo "Wave0 is interactive (env + sync). Use scripts/lambda_bootstrap.sh"
    ;;
  m1)
    run_bg m1_alphamissense "$PY" scripts/join_alphamissense.py \
      --variant-list data/modeling/_tmp_ldblock_associated_sample.parquet \
      --out data/annotations/alphamissense_grch38.parquet \
      --min-free-gb 0
    ;;
  m2_download)
    run_bg m2_vcf_download "$PY" scripts/download_1000g_grch38_all_chrs.py --chroms 1-7
    ;;
  m2_rebuild)
    run_bg m2_rebuild "$PY" scripts/rebuild_chr1_7_score_pgens.py \
      --chroms 1,2,3,4,5,6,7 \
      --memory-mb "$BTA_PLINK_MEMORY_MB" \
      --threads "$BTA_PLINK_THREADS" \
      --jobs "$BTA_REBUILD_JOBS"
    ;;
  m2_downstream)
    run_bg m2_downstream "$PY" scripts/run_genomewide_downstream.py \
      --step score,atlas,intervention \
      --tag genomewide --chroms 1-22 --score-chroms 1-22
    ;;
  gpu)
    run_bg gpu_shap "$PY" scripts/run_shap_attribution.py --tag genomewide || true
    ;;
  help|*)
    echo "Usage: $0 {wave0|m1|m2_download|m2_rebuild|m2_downstream|gpu}"
    exit 1
    ;;
esac
