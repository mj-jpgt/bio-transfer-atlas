# Methods

## Study design

The Biological Transferability Atlas (BTA) predicts cross-ancestry transferability failure
at the **variant** level and evaluates **score-edit interventions** on Catalog PGS scored
in 1000 Genomes Project GRCh38 high-coverage samples (n≈3202). Discovery GWAS labels and
effect sizes primarily use Pan-UK Biobank ancestry-stratified summary statistics for
type 2 diabetes (T2D), coronary artery disease (CAD), body-mass index (BMI), and LDL
cholesterol, with Catalog PGS for the same traits plus white-blood-cell count (WBC),
rheumatoid arthritis (RA), and inflammatory bowel disease (IBD) for positive / stress
controls.

Analyses are split into **variant-scale** (predict high-\(I^2\) per variant) and
**trait-scale** (predict aggregate MAD / high-\(I^2\) rate from trait×pair summaries).
Trait-constant concordance features are **not** peers of AF/LD in variant ablation.

## Genotypes and PGS scoring

Autosomal genotypes were converted to PLINK2 pfiles (`chr{1–22}.score`) and used for
genome-wide scoring of harmonized GRCh38 Catalog scores. Harmonization matches alleles
to the 1000G panel; unmatched variants are dropped. Score matrices are stored as
`score_matrix_grch38_genomewide_genomewide.parquet`.

## Portability labels and features

Variants associated in meta-analysis form the modeling universe. The primary failure label
is high heterogeneity (\(I^2\)) across ancestries (`y_high_I2`). Feature groups include:

- **AF_LD_SEL**: allele-frequency differences, LD scores / tagging metrics, selection proxies
  (nested AF / LD / SEL ablations reported separately)
- **FST**: coarse population differentiation
- **POP_DISTANCE**: ancestry-distance features
- **TRAIT_CONSTANT_Z_DIAGNOSTIC** (not a peer): trait-level Z-score concordance joined as a
  constant within trait—diagnostic only for between-trait rate differences
- Optional: VEP / AlphaMissense / subpopulation AF (within-AFR / within-EUR)

## Predictive models and splits

Primary classifier: histogram gradient boosting (HGB), with XGBoost (CUDA) and a multilayer
perceptron as GPU companions. Evaluation uses **LD-block holdout** (primary claim) and
variant holdout (upper bound; can leak via LD). Metrics: AUROC, AUPRC, calibration.
Paired ΔAUROC vs FST / POP_DISTANCE uses LD-block bootstrap CIs.

## Cross-ancestry Z-score concordance (not Popcorn \(r_g\))

For each trait ∈ {T2D, CAD, BMI, LDL} and ancestry pairs EUR–AFR and EUR–EAS we:

1. Assemble Pan-UKB ancestry Z-statistics across available chromosomes
2. Prefer Popcorn `compute` + `fit` (`--gen_effect`) when inputs succeed; otherwise report
   **Pearson Z-score concordance** with method tag `panukbb_z_concordance_*`
3. Store `z_concordance_by_trait.parquet` (legacy alias `rg_real_by_trait.parquet`) and
   `results/tables/z_concordance_by_trait_pair.csv` / `popcorn_rg_summary.csv`

Trait-scale predictors of mean MAD / high-\(I^2\) rate use concordance, FST summaries,
n variants, and trait class (`trait_scale_portability.csv`).

## Fine-mapping (PolyFun / SuSiE)

Per LD block, `susieR::susie_rss` uses a **signed** 1000G LD matrix from
`plink2 --r-unphased square` (`fallback_mode=signed_ld`), with symmetrization, ridge, and
PSD repair. Residual variance is not estimated under reference LD.
Unsigned |r| from `--r2-unphased`, identity \(R\), and |z|-weights are
`unsigned_abs_r` / `identity_ld` / `absz_weight` and **excluded from primary
`fine_mapped`**. Primary tiers require `fallback_mode == signed_ld` and `in_cs`. Successful
signed extracts record `ld_type=signed_r`, `plink_flag=--r-unphased square`,
`min_eigenvalue_before`, `ridge_added`, and `variant_order_verified` (negative-entry
heuristic only; not a hard assert).

## LD-graph GAT

Within each LD block, variants are nodes; edges connect variants within 250 kb (or
feature-space kNN fallback). Node features = AF_LD_SEL. A 2-layer GAT
(PyTorch Geometric) is trained under LD-block CV and compared to HGB/XGB on the same split.
If GAT ≤ trees, trees remain the primary classifier. Attention weights are **descriptive
graph visualizations**, not causal mechanism evidence; prefer SHAP / grouped permutation
(`grouped_permutation_importance_ldblock.csv`) for attribution claims.

## Interventions (ancestry mean separation)

Modes: risk filters (top 5/10/20%), linear/exponential reweighting, \(F_{ST}\)/MAF/LD
pruning, random drop (negative control), and **duffy_gate** (drop ACKR1 ±200 kb around
rs2814778). Edited weights are re-scored genome-wide.

Primary outcome is **mean absolute EUR–non-EUR score separation (MAD)**—not phenotype
accuracy. We also report variant retention, |β| mass retained, within-ancestry score
variance retained, corr(edited, original), matched-n / matched-mass Monte Carlo (≥500
draws; empirical \(p\) in `intervention_matched_random_controls.csv`), leave-one-score-out
MAD by mode (`intervention_loso_mad_by_mode.csv`), and a retention–MAD curve.

## Concordance analyses (internal vs external)

**Internal Pan-UKB sensitivity:** EUR–AFR β concordance before/after dropping high |Δβ|
variants, with matched random drop of the same *n*. Expected improvement by construction;
labeled `internal_only`—not external validation.

**External PAGE:** GRCh38 liftover of PAGE LDL associations; allele QC ladder
(match/swap/complement; drop ambiguous strand and indels); filter frozen on Pan-UKB only
with assert that PAGE columns are disjoint from filter-training features; require
`n_variants≥500` for `ok_external`; emit `external_page_qc_counts.csv` plus Pearson/
Spearman, sign concordance, slope, and standardized β RMSE. Held-out-trait ranking:
`external_heldout_trait_ranking.csv`.

## Duffy / autoimmune controls

WBC PGS (PGS000191) is stratified by **Duffy-null allele (C) dosage** at rs2814778 after
explicit PLINK allele-orientation audit (`duffy_allele_audit.csv`). Primary AFR contrast:
**dose==2 vs dose&lt;2** with bootstrap CIs; 0/1/2 strata keep n&lt;10 underpowered flags.
ACKR1 score decomposition (lead SNP / window rest / outside) reports `ackr1_fraction`
(`duffy_ackr1_score_decomposition.csv`). RA (PGS004133) and IBD (PGS001288) probe
MHC-heavy scores.

## Software

Python 3 (pandas, scikit-learn, xgboost, torch, torch_geometric), PLINK2, R/`susieR`,
optional Popcorn/PolyFun clones under `tools/`. Reproducible Lambda launchers live in
`scripts/lambda_*.sh`.
