#!/usr/bin/env bash
cd /lambda/nfs/geeg/fairness
source .venv/bin/activate
source scripts/lambda_env.sh
python - <<'PY'
import sys
sys.path.insert(0, "scripts")
from intervention_common import PLINK2
print("PLINK2=", PLINK2)
PY
nohup python scripts/run_genomewide_downstream.py --step score,atlas,intervention --tag genomewide --memory-mb 24000 > results/logs/lambda/m2_downstream.log 2>&1 &
echo DOWN=$!
sleep 20
head -50 results/logs/lambda/m2_downstream.log
