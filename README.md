# Biological Transferability Atlas

Predict **cross-ancestry polygenic score / GWAS portability failure** from open allele-frequency and LD features, then test **score-edit interventions** on 1000 Genomes (GRCh38).

Labels come from multi-ancestry summary-statistic concordance (high \(I^2\), sign discordance)—no individual phenotypes required for the primary analyses.

## Results

Under **LD-block** cross-validation (~109k test variants):

| Analysis | Result |
|----------|--------|
| AF + LD + selection features | AUROC **0.627** for high-\(I^2\) failure |
| AF alone → AF+LD | **0.619 → 0.629** |
| vs coarse \(F_{ST}\) | Higher (paired ΔAUROC CI excludes 0) |
| Feature family permutation | Largest AUROC drop from **AF**, then **LD** |
| Sign-discordance endpoint | AUROC **0.77** |
| Score edits (Catalog PGS) | Frequency-aware filters (\(F_{ST}\)/MAF) reduce ancestry mean score separation (MAD); matched-mass Monte Carlo supports several modes |
| Duffy / WBC positive control | AFR null-homozygous vs dose&lt;2: Δmean PGS **0.119** [0.024, 0.206] |
| PAGE LDL (external) | **4,384** variants after allele QC; β correlation ≈ **0.01** |

Full CSV tables: [`bio-transfer-atlas/results/tables/`](bio-transfer-atlas/results/tables/) · pack: [`bio-transfer-atlas/BTA_robustness_freeze_results_20260720/`](bio-transfer-atlas/BTA_robustness_freeze_results_20260720/)  
Write-up: [`bio-transfer-atlas/paper/`](bio-transfer-atlas/paper/) · [`bio-transfer-atlas/results/BTA_M4_M6_RESULTS_REPORT.md`](bio-transfer-atlas/results/BTA_M4_M6_RESULTS_REPORT.md)

## Setup

```bash
cd bio-transfer-atlas
mamba env create -f environment.yml
mamba activate bta
make test
```

Genotypes, Pan-UKB dumps, and other large inputs are **not** in git. Place them under `bio-transfer-atlas/data/` following the manifests in that tree (or use your existing processed mirrors) before re-running scoring / Duffy scripts.

## Reproduce key analyses

All commands below are from `bio-transfer-atlas/` with the `bta` env active. They write CSVs under `results/tables/`.

```bash
# Nested AF / LD / SEL ablation + paired ΔAUROC
python scripts/run_nested_ablation_and_paired_auroc.py

# Grouped permutation importance (AF / LD / SEL / distance / FST)
python scripts/run_grouped_permutation_importance.py

# Sign-discordance + power sensitivity tables
python scripts/eval_sign_discordance_and_power.py

# Intervention MAD evaluation (uses pre-scored matrices when present)
python scripts/evaluate_intervention.py

# Matched-n / matched-mass Monte Carlo (500 draws) + leave-one-score-out
python scripts/run_intervention_matched_random_controls.py --n-draws 500

# Duffy-null genotype strata + ACKR1 score decomposition (needs chr1 score pfiles)
python scripts/run_duffy_positive_control.py

# PAGE external concordance (uses cached join if present; else GWAS Catalog API)
python scripts/run_external_page_grch38.py --from-joined data/raw/external_sumstats/page_panukbb_joined_grch38.parquet
```

Optional chr22 smoke pipeline:

```bash
make smoke
```

## Repository layout

```
bio-transfer-atlas/
├── configs/          # experiment YAML
├── scripts/          # analysis entrypoints
├── src/bta/          # library code
├── workflow/         # Snakemake rules
├── paper/            # manuscript
├── results/tables/   # result CSVs
├── tests/
└── environment.yml
```

## Data sources

| Source | Use |
|--------|-----|
| 1000 Genomes / IGSR (GRCh38) | Genotypes, ancestry labels |
| PGS Catalog | Score weights |
| Pan-UK Biobank | Multi-ancestry sumstats → labels |
| PAGE (GWAS Catalog) | External β concordance |

## License

Code: [MIT](LICENSE). Upstream genomic datasets keep their own terms.
