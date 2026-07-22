#!/usr/bin/env bash
# Parallel CPU scoring (resume via cached .sscore) + GPU lane on A100.
set -euo pipefail
cd /lambda/nfs/geeg/fairness
source .venv/bin/activate
source scripts/lambda_env.sh
export BTA_SCORE_JOBS="${BTA_SCORE_JOBS:-4}"
export BTA_PLINK_THREADS="${BTA_PLINK_THREADS:-6}"
export BTA_PLINK_MEMORY_MB="${BTA_PLINK_MEMORY_MB:-12000}"
mkdir -p results/logs/lambda

# Stop serial scorer; restart parallel (skips existing .sscore)
pkill -f 'scripts/score_genomewide.py' || true
pkill -f 'run_genomewide_downstream.py' || true
sleep 2

nohup python scripts/run_gpu_lane.py --max-rows 400000 --device cuda \
  > results/logs/lambda/gpu_lane.log 2>&1 &
echo GPU_LANE_PID=$!

nohup bash -c '
  set -e
  python scripts/score_genomewide.py \
    --chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
    --memory-mb '"$BTA_PLINK_MEMORY_MB"' \
    --threads '"$BTA_PLINK_THREADS"' \
    --jobs '"$BTA_SCORE_JOBS"' \
    --out data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet
  python scripts/run_genomewide_downstream.py --step atlas,intervention --tag genomewide --memory-mb 24000
' > results/logs/lambda/m2_downstream.log 2>&1 &
echo SCORE_PID=$!

sleep 8
echo '--- nvidia ---'
nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv
echo '--- procs ---'
pgrep -af 'run_gpu_lane|score_genomewide|plink2' | head -20
echo '--- gpu log ---'
head -30 results/logs/lambda/gpu_lane.log || true
echo '--- score log ---'
tail -20 results/logs/lambda/m2_downstream.log || true
