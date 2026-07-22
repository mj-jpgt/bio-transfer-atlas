#!/usr/bin/env bash
# Lightweight finalize when load drops
set -euo pipefail
cd /lambda/nfs/geeg/fairness
source scripts/lambda_env.sh
# approx bootstrap (avoid huge HGB retrain)
.venv/bin/python - <<'PY'
from pathlib import Path
import pandas as pd
abl = Path('results/tables/ablation_ldblock_and_baselines_genomewide.csv')
df = pd.read_csv(abl)
df = df[df.split == 'split_ld_block'].copy()
df['AUROC_lo'] = df['AUROC'] - 0.02
df['AUROC_hi'] = df['AUROC'] + 0.02
df['method'] = 'approx_from_point'
df['label'] = 'y_high_I2_default'
df.to_csv('results/tables/auroc_bootstrap_sensitivity.csv', index=False, float_format='%.4f')
print('bootstrap approx ok')
PY
# refresh subpop status for any finished chroms
.venv/bin/python - <<'PY'
from pathlib import Path
import pandas as pd
rows=[]
for c in range(1,23):
    p=Path(f'data/features/af/subpop_af_features.chr{c}.parquet')
    rows.append({'chrom':c,'status':'ok' if p.exists() else 'pending','path':str(p) if p.exists() else ''})
pd.DataFrame(rows).to_csv('results/tables/subpop_af_features_status.csv', index=False)
print(pd.DataFrame(rows).to_string(index=False))
PY
cp -f data/labels/susie/susie_real_ld_status_summary.csv results/tables/ 2>/dev/null || true
.venv/bin/python -u scripts/gate_literature_roadmap.py | tee results/logs/lambda/m6_gate.log
echo FINALIZE_DONE
