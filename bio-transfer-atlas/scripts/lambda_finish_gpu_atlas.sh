#!/usr/bin/env bash
# Finish GPU torch MLP + genomewide atlas/intervention on Lambda A100.
set -euo pipefail
ROOT=/lambda/nfs/geeg/fairness
cd "$ROOT"
# shellcheck disable=SC1091
source scripts/lambda_env.sh
mkdir -p results/logs/lambda

export OMP_NUM_THREADS="${BTA_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${BTA_THREADS:-8}"
export MKL_NUM_THREADS="${BTA_THREADS:-8}"
export NUMEXPR_NUM_THREADS="${BTA_THREADS:-8}"

echo "[finish] torch check"
.venv/bin/python - <<'PY'
import torch
print(torch.__version__, torch.cuda.is_available(), torch.cuda.get_device_name(0))
print("VRAM", round(torch.cuda.get_device_properties(0).total_memory/1024**3,1), "GB")
PY

echo "[finish] launch atlas+intervention (CPU/RAM heavy)"
nohup .venv/bin/python -u scripts/run_genomewide_downstream.py \
  --step atlas,intervention \
  --tag genomewide \
  --memory-mb "${BTA_MEMORY_MB:-24000}" \
  > results/logs/lambda/m2_atlas_intervention.log 2>&1 &
echo "ATLAS_PID=$!"

echo "[finish] launch sequential GPU lane: torch MLP then heavy XGB (avoid VRAM contention)"
nohup bash -c '
  set -e
  cd /lambda/nfs/geeg/fairness
  source scripts/lambda_env.sh
  .venv/bin/python -u scripts/run_gpu_lane.py \
    --only torch \
    --max-rows 800000 \
    --torch-hidden 4096 \
    --torch-epochs 25 \
    --device cuda \
    > results/logs/lambda/gpu_torch_mlp.log 2>&1
  .venv/bin/python -u scripts/run_gpu_lane.py \
    --only xgb \
    --max-rows 600000 \
    --device cuda \
    > results/logs/lambda/gpu_xgb_heavy.log 2>&1
  echo GPU_SEQ_DONE >> results/logs/lambda/gpu_torch_mlp.log
' >/dev/null 2>&1 &
echo "GPU_SEQ_PID=$!"

sleep 15
nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv
pgrep -af 'run_genomewide_downstream|run_gpu_lane' | grep -v grep || true
echo "[finish] launched; tails:"
tail -n 8 results/logs/lambda/gpu_torch_mlp.log || true
tail -n 8 results/logs/lambda/m2_atlas_intervention.log || true
