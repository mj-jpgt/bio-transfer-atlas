#!/usr/bin/env bash
cd /lambda/nfs/geeg/fairness
source scripts/lambda_env.sh
.venv/bin/pip install -q torch_geometric
nohup .venv/bin/python -u scripts/train_ld_gat.py --device cuda --epochs 25 --max-blocks 400 \
  > results/logs/lambda/m5_gat.log 2>&1 &
echo GAT_PID=$!
.venv/bin/python -u scripts/run_external_validation_ci.py \
  > results/logs/lambda/m4_external_ci.log 2>&1
echo EXTERNAL_CI_DONE
