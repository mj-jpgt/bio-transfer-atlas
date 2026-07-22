#!/usr/bin/env bash
# One-shot bootstrap on Lambda: venv, plink2 linux, R/susieR, ldsc, Popcorn, polyfun.
set -euo pipefail
ROOT="${BTA_ROOT:-/lambda/nfs/geeg/fairness}"
cd "$ROOT"
# shellcheck disable=SC1091
source "$ROOT/scripts/lambda_env.sh" 2>/dev/null || true
mkdir -p "$ROOT/tools" "$ROOT/data/raw/ldscores" "$ROOT/results/logs/lambda" "$TMPDIR"

echo "=== Python venv ==="
if [[ ! -d .venv ]]; then
  python3 -m venv .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate
pip install -U pip wheel
pip install numpy pandas pyarrow scipy scikit-learn statsmodels matplotlib seaborn \
  requests tqdm pyyaml joblib shap xgboost lightgbm pandera pydantic typer rich loguru \
  gseapy 2>&1 | tail -20

echo "=== Linux plink2 ==="
if [[ ! -x tools/plink2/plink2 ]]; then
  mkdir -p tools/plink2
  cd /tmp
  curl -fLsS -o plink2.zip "https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_latest.zip" \
    || curl -fLsS -o plink2.zip "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20241114.zip"
  unzip -o plink2.zip -d "$ROOT/tools/plink2"
  chmod +x "$ROOT/tools/plink2/plink2"
  cd "$ROOT"
fi
"$ROOT/tools/plink2/plink2" --version | head -2

echo "=== R + susieR ==="
if ! command -v Rscript >/dev/null 2>&1; then
  sudo DEBIAN_FRONTEND=noninteractive apt-get update -qq
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev
fi
Rscript -e 'if (!requireNamespace("susieR", quietly=TRUE)) install.packages(c("susieR","Rcpp","RcppArmadillo"), repos="https://cloud.r-project.org")'
Rscript -e 'cat("susieR=", requireNamespace("susieR", quietly=TRUE), "\n")'

echo "=== ldsc ==="
if [[ ! -d tools/ldsc ]]; then
  git clone --depth 1 https://github.com/bulik/ldsc.git tools/ldsc
fi
pip install -q bitarray 2>/dev/null || true

echo "=== Popcorn ==="
if [[ ! -d tools/Popcorn ]]; then
  git clone --depth 1 https://github.com/brielin/Popcorn.git tools/Popcorn
fi

echo "=== PolyFun ==="
if [[ ! -d tools/polyfun ]]; then
  git clone --depth 1 https://github.com/omerwe/polyfun.git tools/polyfun
fi

echo "=== LD score tarballs (EUR baseline; AFR/EAS if available) ==="
LDS="$ROOT/data/raw/ldscores"
mkdir -p "$LDS"
cd "$LDS"
if [[ ! -f eur_w_ld_chr.tar.bz2 && ! -d eur_w_ld_chr ]]; then
  curl -fLsS -O https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 || true
  tar -xjf eur_w_ld_chr.tar.bz2 2>/dev/null || true
fi
cd "$ROOT"

echo "BOOTSTRAP_OK"
plink2 --version | head -1
python -c "import pandas,sklearn,shap; print('py_ok', pandas.__version__)"
Rscript -e 'cat(as.character(packageVersion("susieR")), "\n")'
