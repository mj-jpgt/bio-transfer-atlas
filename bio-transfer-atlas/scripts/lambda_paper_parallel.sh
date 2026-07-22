#!/usr/bin/env bash
# Parallel paper campaign: long lanes in background; M4 scoring first in foreground path.
set -euo pipefail
ROOT=/lambda/nfs/geeg/fairness
cd "$ROOT"
source scripts/lambda_env.sh
mkdir -p results/logs/lambda results/tables results/figures paper
PY=.venv/bin/python
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
log() { echo "[$(date -Is)] $*"; }

log "Download WBC/RA/IBD PGS"
$PY - <<'PY'
from pathlib import Path
import requests
ids = ["PGS000191", "PGS004133", "PGS001288"]
base = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs}/ScoringFiles/Harmonized/{pgs}_hmPOS_GRCh38.txt.gz"
for pgs in ids:
    dest = Path(f"data/raw/pgs_catalog/scores/{pgs}/{pgs}_hmPOS_GRCh38.txt.gz")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(pgs, "exists", dest.stat().st_size); continue
    r = requests.get(base.format(pgs=pgs), timeout=180); r.raise_for_status()
    dest.write_bytes(r.content); print(pgs, "ok", len(r.content))
PY

log "Harmonize"
$PY -u scripts/harmonize_pgs_genomewide.py > results/logs/lambda/m4_harmonize.log 2>&1

log "Score expansion PGS"
nohup $PY -u scripts/score_genomewide.py \
  --pgs-ids PGS000191,PGS004133,PGS001288 \
  --chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
  --jobs 8 --threads 3 --memory-mb 6000 \
  --out data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet \
  > results/logs/lambda/m4_score_expand.log 2>&1 &
echo SCORE_PID=$!

log "Launch subpop chr1-7"
nohup $PY -u scripts/compute_subpop_af_features.py \
  --chroms 1,2,3,4,5,6,7 --jobs 4 --memory-mb 8192 \
  > results/logs/lambda/m5_subpop_chr1_7.log 2>&1 &
echo SUBPOP_PID=$!

log "Launch popcorn autosome"
nohup $PY -u scripts/run_popcorn_rg.py --autosome \
  > results/logs/lambda/m5_popcorn_autosome.log 2>&1 &
echo POPCORN_PID=$!

log "Launch PAGE external"
nohup $PY -u scripts/run_external_page_lean.py --chrom 22 --size 2000 \
  > results/logs/lambda/m4_external_page.log 2>&1 &
echo PAGE_PID=$!

log "Launch SuSiE real-LD (priority chroms + traits)"
nohup bash -c '
  source scripts/lambda_env.sh
  export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
  for trait in T2D CAD BMI LDL; do
    for anc in EUR AFR; do
      for chrom in 1 2 3 6 7 19 22; do
        .venv/bin/python -u scripts/run_polyfun_susie.py \
          --trait "$trait" --chrom "$chrom" --anc "$anc" \
          --jobs 6 --max-blocks 30 --prefer-real-ld \
          >> results/logs/lambda/m5_susie_farm.log 2>&1 || true
      done
    done
  done
  .venv/bin/python -u scripts/_postprocess_susie_cs.py >> results/logs/lambda/m5_susie_farm.log 2>&1 || true
  .venv/bin/python -u scripts/build_finemap_tier_labels.py --tag genomewide_susie >> results/logs/lambda/m5_susie_farm.log 2>&1 || true
  .venv/bin/python -u scripts/eval_finemap_tiers_lean.py >> results/logs/lambda/m5_susie_farm.log 2>&1 || true
' > results/logs/lambda/m5_susie_wrapper.log 2>&1 &
echo SUSIE_PID=$!

log "Wait for scoring to finish before Duffy/intervene/GAT"
wait ${SCORE_PID}

log "Duffy control"
$PY -u scripts/run_duffy_positive_control.py > results/logs/lambda/m4_duffy.log 2>&1 || true

log "Apply interventions"
$PY -u scripts/apply_intervention.py \
  --modes filter_10,random,fst,maf,duffy_gate \
  > results/logs/lambda/m4_apply_int.log 2>&1 || true

log "Score key interventions for new PGS"
nohup $PY -u scripts/score_intervention.py \
  --pgs-ids PGS000191,PGS004133,PGS001288 \
  --modes fst,maf,duffy_gate,random,filter_10 \
  --chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
  --jobs 6 --threads 2 --memory-mb 5000 --tag genomewide \
  > results/logs/lambda/m4_score_int.log 2>&1 &
echo SCORE_INT_PID=$!

log "GAT on A100 (uses existing LD-block sample)"
nohup $PY -u scripts/train_ld_gat.py --device cuda --epochs 25 --max-blocks 400 \
  > results/logs/lambda/m5_gat.log 2>&1 &
echo GAT_PID=$!

log "Wait popcorn then re-ablate"
wait ${POPCORN_PID} || true
$PY -u scripts/run_ldblock_and_baselines.py > results/logs/lambda/m5_reablate.log 2>&1 || true

wait ${SCORE_INT_PID} || true
$PY -u scripts/evaluate_intervention.py --tag genomewide \
  > results/logs/lambda/m4_eval_int.log 2>&1 || true

wait ${GAT_PID} || true
wait ${SUBPOP_PID} || true
wait ${PAGE_PID} || true
wait ${SUSIE_PID} || true

log "Figures + robustness + gate"
$PY -u scripts/figures/make_paper_figures.py > results/logs/lambda/m5_figures.log 2>&1 || true
$PY -u scripts/run_robustness_bootstrap.py > results/logs/lambda/m5_robust.log 2>&1 || true
$PY -u scripts/gate_literature_roadmap.py > results/logs/lambda/m6_gate.log 2>&1 || true

log "PARALLEL_CAMPAIGN_DONE"
cat results/logs/lambda/m6_gate.log || true
