#!/usr/bin/env bash
# Full paper-completion campaign M4→M6 on Lambda A100 host.
# Usage: bash scripts/lambda_paper_m4_m6.sh
set -euo pipefail
ROOT=/lambda/nfs/geeg/fairness
cd "$ROOT"
# shellcheck disable=SC1091
source scripts/lambda_env.sh
mkdir -p results/logs/lambda results/tables results/figures paper
PY=.venv/bin/python
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"

log() { echo "[$(date -Is)] $*"; }

# ── Download new PGS if missing ──────────────────────────────────────────────
log "Download expansion PGS (WBC/RA/IBD)"
$PY - <<'PY'
from pathlib import Path
import requests
ROOT = Path(".")
ids = ["PGS000191", "PGS004133", "PGS001288"]
base = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs}/ScoringFiles/Harmonized/{pgs}_hmPOS_GRCh38.txt.gz"
for pgs in ids:
    dest = ROOT / f"data/raw/pgs_catalog/scores/{pgs}/{pgs}_hmPOS_GRCh38.txt.gz"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(pgs, "exists", dest.stat().st_size)
        continue
    url = base.format(pgs=pgs)
    print("GET", url)
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    dest.write_bytes(r.content)
    print(pgs, "ok", len(r.content))
PY

# ── M4: harmonize + score expansion ─────────────────────────────────────────
log "Harmonize genomewide PGS"
$PY -u scripts/harmonize_pgs_genomewide.py > results/logs/lambda/m4_harmonize.log 2>&1 || {
  echo "harmonize failed; see m4_harmonize.log"; tail -50 results/logs/lambda/m4_harmonize.log; exit 1;
}

log "Score WBC/RA/IBD into genomewide matrix"
$PY -u scripts/score_genomewide.py \
  --pgs-ids PGS000191,PGS004133,PGS001288 \
  --chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
  --jobs 8 --threads 3 --memory-mb 6000 \
  --out data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet \
  > results/logs/lambda/m4_score_expand.log 2>&1

# ── M4a Duffy control (needs WBC scores) ────────────────────────────────────
log "Duffy positive control"
$PY -u scripts/run_duffy_positive_control.py > results/logs/lambda/m4_duffy.log 2>&1 || true

# ── Apply interventions including duffy_gate for all PGS ────────────────────
log "Apply interventions (incl duffy_gate)"
$PY -u scripts/apply_intervention.py \
  --modes filter_5,filter_10,filter_20,reweight_linear,reweight_exp,flag,random,fst,ld,maf,duffy_gate \
  > results/logs/lambda/m4_apply_int.log 2>&1

log "Score interventions for expansion PGS (parallel)"
$PY -u scripts/score_intervention.py \
  --pgs-ids PGS000191,PGS004133,PGS001288 \
  --modes fst,maf,duffy_gate,random,filter_10 \
  --chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
  --jobs 6 --threads 2 --memory-mb 5000 \
  --tag genomewide \
  > results/logs/lambda/m4_score_int.log 2>&1

log "Evaluate interventions"
$PY -u scripts/evaluate_intervention.py --tag genomewide \
  > results/logs/lambda/m4_eval_int.log 2>&1 || true

# ── M4c external PAGE ───────────────────────────────────────────────────────
log "PAGE external validation"
$PY -u scripts/run_external_page_lean.py --chrom 22 --size 2000 \
  > results/logs/lambda/m4_external_page.log 2>&1 || true

# ── M5 subpop chr1-7 (background-friendly but wait here for completeness) ───
log "Subpop AF chr1-7"
$PY -u scripts/compute_subpop_af_features.py \
  --chroms 1,2,3,4,5,6,7 --jobs 4 --memory-mb 8192 \
  > results/logs/lambda/m5_subpop_chr1_7.log 2>&1 || true

# ── M5 Popcorn autosomal ────────────────────────────────────────────────────
log "Popcorn / Z-corr autosome"
$PY -u scripts/run_popcorn_rg.py --autosome --skip-download \
  > results/logs/lambda/m5_popcorn_autosome.log 2>&1 || \
$PY -u scripts/run_popcorn_rg.py --autosome \
  > results/logs/lambda/m5_popcorn_autosome.log 2>&1 || true

log "Re-ablate LD-block with refreshed RG_REAL"
$PY -u scripts/run_ldblock_and_baselines.py \
  > results/logs/lambda/m5_reablate.log 2>&1 || true

# ── M5 SuSiE real-LD farm (traits x EUR/AFR, chroms 1-22 capped) ────────────
log "SuSiE real-LD farm"
for trait in T2D CAD BMI LDL; do
  for anc in EUR AFR; do
    for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
      $PY -u scripts/run_polyfun_susie.py \
        --trait "$trait" --chrom "$chrom" --anc "$anc" \
        --jobs 8 --max-blocks 40 --prefer-real-ld \
        >> results/logs/lambda/m5_susie_farm.log 2>&1 || true
    done
  done
done
$PY -u scripts/_postprocess_susie_cs.py >> results/logs/lambda/m5_susie_farm.log 2>&1 || true
$PY -u scripts/build_finemap_tier_labels.py --tag genomewide_susie >> results/logs/lambda/m5_susie_farm.log 2>&1 || true
$PY -u scripts/eval_finemap_tiers_lean.py >> results/logs/lambda/m5_susie_farm.log 2>&1 || true

# ── M5 GAT on A100 ──────────────────────────────────────────────────────────
log "Train LD-block GAT"
$PY -u scripts/train_ld_gat.py --device cuda --epochs 25 --max-blocks 400 \
  > results/logs/lambda/m5_gat.log 2>&1 || true

# ── M5 figures + robustness ─────────────────────────────────────────────────
log "Figures + robustness"
$PY -u scripts/figures/make_paper_figures.py > results/logs/lambda/m5_figures.log 2>&1 || true
$PY -u scripts/run_robustness_bootstrap.py > results/logs/lambda/m5_robust.log 2>&1 || true

# ── M6 paper gate ───────────────────────────────────────────────────────────
log "Literature roadmap gate"
$PY -u scripts/gate_literature_roadmap.py > results/logs/lambda/m6_gate.log 2>&1 || true

log "CAMPAIGN_DONE"
tail -30 results/logs/lambda/m6_gate.log || true
