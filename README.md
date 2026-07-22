# Biological Transferability Atlas

Predict **cross-ancestry polygenic score / GWAS portability failure** from open allele-frequency and LD features, then test **score-edit interventions** on 1000 Genomes (GRCh38).

Labels come from multi-ancestry summary-statistic concordance (high \(I^2\), sign discordance)—no individual phenotypes required for the primary analyses.

## Abstract

Cross-ancestry polygenic scores often lose accuracy outside European discovery cohorts. We build an open-data **Biological Transferability Atlas**: variant-level models that predict GWAS/PRS portability failure from allele-frequency (AF), LD, and selection features, then test score-edit interventions on 1000 Genomes GRCh38 genotypes. Under LD-block cross-validation, AF+LD features achieve AUROC ≈ **0.627**, exceeding coarse \(F_{ST}\); population-distance features are competitive (~0.617) and should not be dismissed as uniformly weak. Frequency-aware interventions (\(F_{ST}\)/MAF filters) reduce mean absolute EUR–non-EUR score gaps on average, whereas pruning top predicted-risk variants often worsens gaps. Positive controls (Duffy/WBC; autoimmune MHC stress-tests) and open external sumstat concordance with bootstrap CIs anchor claim discipline. Fine-mapping tiers and an LD-graph GAT provide descriptive attention weights without claiming that SuSiE or graph nets “solve” allele-frequency bias.

## Contributions

1. **Variant-level portability-risk models** under LD-block cross-validation, using AF, LD, and selection features to predict high-\(I^2\) GWAS discordance (and related endpoints such as sign discordance).
2. **Honest peer comparisons**: AF/LD vs coarse \(F_{ST}\) and population-distance encodings at the variant scale; trait-level Z-score concordance kept as a trait-scale analysis, not an AF/LD peer contest.
3. **Score-edit intervention bake-off** on Catalog PGS applied to 1000 Genomes GRCh38, scored as ancestry **mean separation** (MAD)—not phenotype accuracy—with matched-mass Monte Carlo and leave-one-score-out controls.
4. **Biological anchors**: Duffy/WBC positive control (allele-audited) and MHC-heavy autoimmune stress-tests; internal Pan-UKB concordance sensitivity vs external PAGE after GRCh38 liftover.
5. **Claim discipline**: SuSiE primary = signed LD only; GAT attention is descriptive graph weighting, not mechanism evidence.

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

Full CSV tables and figures: [`bio-transfer-atlas/BTA_robustness_freeze_results_20260720/`](bio-transfer-atlas/BTA_robustness_freeze_results_20260720/)  
Narrative report: [`bio-transfer-atlas/BTA_robustness_freeze_results_20260720/reports/BTA_M4_M6_RESULTS_REPORT.md`](bio-transfer-atlas/BTA_robustness_freeze_results_20260720/reports/BTA_M4_M6_RESULTS_REPORT.md)

## Setup

```bash
cd bio-transfer-atlas
mamba env create -f environment.yml
mamba activate bta
make test
```

Genotypes, Pan-UKB dumps, and other large inputs are **not** in git. Place them under `bio-transfer-atlas/data/` following the manifests in that tree (or use your existing processed mirrors) before re-running scoring / Duffy scripts.

## Reproduce key analyses

From `bio-transfer-atlas/` with the `bta` env active. Scripts write local outputs under `results/tables/`; published tables are in `BTA_robustness_freeze_results_20260720/tables/`.

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

## Repository layout

```
bio-transfer-atlas/
├── configs/          # experiment YAML
├── scripts/          # analysis entrypoints
├── src/bta/          # library code
├── BTA_robustness_freeze_results_20260720/  # result tables + figures
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
