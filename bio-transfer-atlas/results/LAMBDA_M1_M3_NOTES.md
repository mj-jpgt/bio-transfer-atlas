# Lambda M1–M3 upgrade notes (2026-07-20)

Work ran on `/lambda/nfs/geeg/fairness` (A100 40GB, ~216 GB RAM).

## Milestone 1 — AlphaMissense / VEP×AF (done)
- Joined GCS `AlphaMissense_hg38.tsv.gz` → `data/annotations/alphamissense_grch38.parquet`
- Re-ran `eval_vep_af_interaction.py` with `vep_source=alphamissense`
- Result: AF_LD_SEL AUROC ≈ 0.730; AF_LD_SEL+VEP_AF ≈ 0.730 (sparse missense overlap)

## Milestone 2 — chr1–7 + genomewide score/atlas/intervention (done)
- Parallel IGSR GRCh38 VCF download + PLINK rebuild
- Score matrix: **3202 × 9 PGS** → `score_matrix_grch38_genomewide_genomewide.parquet`
- Portability model: Test AUROC (split_variant) **0.6869**
- Intervention apply + **parallel PLINK rescoring** (12 jobs × 3 threads; 1980 tasks, resume-friendly)
- Eval: `intervention_results.genomewide.csv` / `intervention_summary.genomewide.txt`
  - Mixed MAD reductions; **fst** often strongest (e.g. CAD PGS large positive reduction); filters sometimes worsen MAD

## Milestone 3A — RG_REAL (done)
- Pan-UKB chr22 ancestry betas → cross-ancestry Z correlation
- LD-block CV: **AF_LD_SEL 0.627 > RG_REAL 0.566** (also > POP_DISTANCE 0.617, RG_PROXY 0.618)
- GPU XGBoost companion: same ordering

## Milestone 3B — SuSiE tiers (done)
- `susieR` + farm T2D/CAD/BMI/LDL × EUR/AFR chr22 (`mode=susie_cs`)
- |z|-weight PIP fallback; CS = top-decile PIP
- `data/labels/susie/susie_*.parquet`, `finemap_tiers_genomewide.parquet`

## Milestone 3C — Subpop AF (done chr8–22)
- Fixed plink2 `.afreq` `#CHROM` header parsing
- All chroms 8–22 parquet written; status → `subpop_af_features_status.csv`
- Within-AFR AF maxdiff means ≈ 0.023–0.024

## GPU lane (done; Torch + expanded XGB)
- Torch 2.6.0+cu124 on A100; sequential torch MLP then heavy XGB (avoid VRAM contention)
- LD-block AF_LD_SEL: **XGB CUDA 0.620** | Torch MLP 0.611 | XGB depth12/1500 0.600 | RFF-MLP 0.549 (overfit)
- Split-variant AF_LD_SEL XGB CUDA **0.826** (leakier split; use LD-block for claims)
- SHAP/gain attribution: dominant **AF** features (`shap_xgboost_gpu_attribution.csv`)
- Tabular models only need ~2–4 GB VRAM; PLINK scoring is CPU/RAM (parallelized instead)
