#!/usr/bin/env bash
set -euo pipefail
cd /lambda/nfs/geeg/fairness
source .venv/bin/activate
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
mkdir -p results/logs/lambda

# Free RAM: stop SHAP if artifact already present
if [[ -f results/tables/shap_mechanism_attribution_genomewide.csv ]]; then
  pkill -f run_shap_attribution.py || true
fi

LOG=results/logs/lambda/m3b_susie_refarm.log
: > "$LOG"
for trait in T2D CAD BMI LDL; do
  for anc in EUR AFR; do
    echo "=== $trait $anc ===" | tee -a "$LOG"
    python scripts/run_polyfun_susie.py --trait "$trait" --chrom 22 --anc "$anc" --jobs 16 --max-blocks 80 2>&1 | tee -a "$LOG"
  done
done

python scripts/build_finemap_tier_labels.py --tag genomewide 2>&1 | tee -a "$LOG"
python scripts/eval_finemap_tiers_lean.py 2>&1 | tee -a "$LOG" || true

python - <<'PY'
import pandas as pd
from pathlib import Path
s = pd.read_parquet("data/labels/susie/susie_T2D_EUR.parquet")
print("T2D_EUR in_cs", int(s.in_cs.sum()), "max_pip", float(s.pip.max()), "n", len(s))
t = pd.read_parquet("data/labels/finemap_tiers_genomewide.parquet")
print(t.finemap_tier.value_counts())
print(t.tier_method.value_counts())
root = Path("data/interim/1000g_grch38")
rows = []
for c in range(1, 8):
    p = root / f"chr{c}.score.pgen"
    rows.append({
        "chrom": c,
        "status": "ok" if p.exists() and p.stat().st_size > 1000 else "pending",
        "bytes": p.stat().st_size if p.exists() else 0,
    })
pd.DataFrame(rows).to_csv("results/tables/score_pgen_chr1_7_rebuild_status.csv", index=False)
print(pd.DataFrame(rows))
PY

# If ablation still running leave it; else print
if pgrep -f run_ldblock_and_baselines.py >/dev/null; then
  echo ABLATION_RUNNING
else
  echo ABLATION_DONE
  tail -40 results/logs/lambda/m3a_ablation.log || true
fi
echo REFARM_DONE
