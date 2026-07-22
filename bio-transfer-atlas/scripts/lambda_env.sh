#!/usr/bin/env bash
# Lambda A100 compute profile + wave launcher for bio-transfer-atlas.
# Usage: source scripts/lambda_env.sh
#        bash scripts/lambda_run_wave.sh wave0|m1|m2|m3|gpu

set -euo pipefail

export BTA_ROOT="${BTA_ROOT:-/lambda/nfs/geeg/fairness}"
cd "$BTA_ROOT"

export BTA_PLINK_MEMORY_MB="${BTA_PLINK_MEMORY_MB:-24000}"
export BTA_PLINK_THREADS="${BTA_PLINK_THREADS:-8}"
export BTA_MIN_FREE_GB="${BTA_MIN_FREE_GB:-0}"
export BTA_REBUILD_JOBS="${BTA_REBUILD_JOBS:-2}"
export BTA_DOWNLOAD_JOBS="${BTA_DOWNLOAD_JOBS:-4}"
export TMPDIR="${TMPDIR:-/tmp/bta}"
mkdir -p "$TMPDIR" "$BTA_ROOT/results/logs/lambda"

# Prefer project venv
if [[ -f "$BTA_ROOT/.venv/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source "$BTA_ROOT/.venv/bin/activate"
fi

# Prefer Linux plink2 in tools/
if [[ -x "$BTA_ROOT/tools/plink2/plink2" ]]; then
  export PATH="$BTA_ROOT/tools/plink2:$PATH"
fi
if [[ -d "$BTA_ROOT/tools/ldsc" ]]; then
  export PATH="$BTA_ROOT/tools/ldsc:$PATH"
  export PYTHONPATH="${PYTHONPATH:-}:$BTA_ROOT/tools/ldsc"
fi
if [[ -d "$BTA_ROOT/tools/Popcorn" ]]; then
  export PYTHONPATH="${PYTHONPATH:-}:$BTA_ROOT/tools/Popcorn"
fi
if [[ -d "$BTA_ROOT/tools/polyfun" ]]; then
  export PYTHONPATH="${PYTHONPATH:-}:$BTA_ROOT/tools/polyfun"
fi

echo "[lambda_env] ROOT=$BTA_ROOT memory=${BTA_PLINK_MEMORY_MB}MB threads=${BTA_PLINK_THREADS}"
