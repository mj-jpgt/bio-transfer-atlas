#!/usr/bin/env bash
# Rebuild score pgens as VCFs finish (pool of 2).
set -euo pipefail
ROOT=/lambda/nfs/geeg/fairness
cd "$ROOT"
# shellcheck disable=SC1091
source scripts/lambda_env.sh
source .venv/bin/activate
VCF_DIR=data/raw/1000g/vcf_grch38
TEMPLATE='1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
LOG=results/logs/lambda/m2_rebuild_watch.log
mkdir -p results/logs/lambda /tmp/bta
export TMPDIR=/tmp/bta

rebuild_one() {
  local chrom=$1
  local score="data/interim/1000g_grch38/chr${chrom}.score.pgen"
  if [[ -f "$score" && $(stat -c%s "$score") -gt 1000 ]]; then
    echo "chr${chrom}: already have score" | tee -a "$LOG"
    return 0
  fi
  echo "chr${chrom}: rebuilding $(date -Is)" | tee -a "$LOG"
  python scripts/rebuild_score_pfile_conservative.py \
    --chrom "$chrom" \
    --memory-mb "${BTA_PLINK_MEMORY_MB:-24000}" \
    --threads "${BTA_PLINK_THREADS:-8}" \
    --min-free-gb 0 \
    --work-root "/tmp/bta/bta_chr${chrom}" \
    --allow-qc-fallback >>"$LOG" 2>&1
  echo "chr${chrom}: exit=$?" | tee -a "$LOG"
}

# Wait until all VCFs present or timeout hours
DEADLINE=$(( $(date +%s) + 3600*12 ))
declare -A STARTED
while true; do
  pending=0
  for chrom in 1 2 3 4 5 6 7; do
    score="data/interim/1000g_grch38/chr${chrom}.score.pgen"
    if [[ -f "$score" && $(stat -c%s "$score") -gt 1000 ]]; then
      continue
    fi
    vcf=$(printf "$TEMPLATE" "$chrom")
    if [[ -f "$VCF_DIR/$vcf" && ! -f "$VCF_DIR/${vcf}.partial" ]]; then
      if [[ -z "${STARTED[$chrom]:-}" ]]; then
        # limit concurrency
        while [[ $(jobs -rp | wc -l) -ge 2 ]]; do sleep 30; done
        STARTED[$chrom]=1
        rebuild_one "$chrom" &
      fi
    else
      pending=$((pending+1))
    fi
  done
  # all done?
  ok=0
  for chrom in 1 2 3 4 5 6 7; do
    score="data/interim/1000g_grch38/chr${chrom}.score.pgen"
    if [[ -f "$score" && $(stat -c%s "$score") -gt 1000 ]]; then
      ok=$((ok+1))
    fi
  done
  if [[ $ok -eq 7 ]]; then
    echo "ALL_CHR1_7_OK" | tee -a "$LOG"
    wait || true
    python scripts/rebuild_chr1_7_score_pgens.py --chroms 1,2,3,4,5,6,7 --memory-mb 24000 --threads 8 --jobs 1 >>"$LOG" 2>&1 || true
    break
  fi
  if [[ $(date +%s) -gt $DEADLINE ]]; then
    echo "TIMEOUT pending=$pending ok=$ok" | tee -a "$LOG"
    wait || true
    break
  fi
  sleep 60
done
echo "WATCH_DONE" | tee -a "$LOG"
