# Biological Transferability Atlas

Polygenic scores (PGS) and GWAS effects often **stop working** when applied outside the ancestries they were trained on. That failure is one of the central equity and translation problems in human genetics: risk tools built mostly in European cohorts can systematically mis-rank people with African, East Asian, South Asian, or admixed ancestry.

This project asks a sharper question than “does accuracy drop?”:

> **Can we predict—*without phenotypes*—which variants and scores are likely to fail across ancestries, what population-genetic mechanism drives that failure, and which score edits actually shrink ancestry gaps?**

The **Biological Transferability Atlas (BTA)** is an open-data answer: a genome-wide map of portability risk built from 1000 Genomes genotypes, PGS Catalog weights, and multi-ancestry summary statistics (Pan-UK Biobank, with PAGE as an external check).

## Why this matters

**For the field.** Much of the literature documents that PGS transfer poorly, or treats ancestry distance as a single dial. Recent work argues that **allele-frequency (AF) and linkage-disequilibrium (LD) structure**—not mysterious effect-size chaos everywhere—drive a large share of the problem. BTA turns that idea into an operational, genome-wide testbed: variant-level failure labels from cross-ancestry GWAS discordance, predictors from AF/LD/selection features, and ablations under **LD-block** cross-validation so claims are not inflated by leakage.

**For translational / clinical genomics.** Before a score is used for screening, triage, or risk communication in a new population, teams need to know (1) which loci are structurally fragile, (2) whether a proposed “fix” merely compresses score distributions or actually targets biology, and (3) when fancy tools (fine-mapping, graph nets) are descriptive lenses rather than cures. BTA stress-tests **score edits** on open genotypes and scores them as ancestry **mean separation**—honest about what we can measure without biobank phenotypes, and useful as a pre-deployment triage layer.

**What we are *not* claiming.** We do not claim improved clinical accuracy, calibrated absolute risk, or that editing a Catalog PGS replaces multi-ancestry discovery. Interventions here reduce score gaps across ancestries; phenotype association remains future work.

## Abstract

Cross-ancestry polygenic scores often lose accuracy outside European discovery cohorts. We build an open-data **Biological Transferability Atlas**: genome-wide models that predict GWAS/PRS portability failure from allele-frequency, LD, and selection features, then test score-edit interventions on 1000 Genomes GRCh38 genotypes. Under LD-block cross-validation (~109k test variants), AF+LD features achieve AUROC ≈ **0.627** for high-\(I^2\) failure—better than coarse \(F_{ST}\), with population-distance features competitive rather than irrelevant. Frequency-aware edits (\(F_{ST}\)/MAF filters) shrink mean EUR–non-EUR score gaps on average, while pruning the highest predicted-risk variants often **widens** them. A Duffy/WBC positive control and MHC-heavy autoimmune scores anchor biology; external PAGE concordance stays near zero after allele QC. Fine-mapping (signed-LD SuSiE) and an LD-graph GAT clarify *where* risk concentrates without claiming to solve AF bias.

## What we built (contributions)

1. **A phenotype-free, genome-wide portability atlas.** End-to-end pipeline on 1000 Genomes GRCh38 + open sumstats: score-shift maps across ancestries, multi-ancestry concordance labels, AF/LD/selection features, and pathway-level risk aggregation—so others can study transfer failure without private biobank access.

2. **Predictive models of *where* transfer fails.** Variant-level models that forecast high cross-ancestry GWAS discordance (and sign discordance) from population-genetic features under rigorous LD-block CV—moving from documenting decay to **locating** fragile loci.

3. **Mechanism-facing evidence, not a single-number ancestry dial.** Nested ablations and grouped permutation show **AF**, then **LD**, carry most of the signal; coarse \(F_{ST}\) underperforms; trait-level Z-score concordance is kept at the **trait scale** (not falsely pitted as a variant-level peer of AF/LD).

4. **Intervention triage that can say “no.”** A bake-off of Catalog PGS edits on real genotypes: which filters reduce ancestry mean separation (MAD), which backfire, with matched-mass Monte Carlo and leave-one-score-out controls—so “drop the risky SNPs” is tested, not assumed.

5. **Biological and external stress tests.** Allele-audited Duffy/WBC positive control (ACKR1), MHC-heavy autoimmune scores, internal Pan-UKB sensitivity checks, and external PAGE after GRCh38 liftover—plus honest limits on SuSiE/GAT as descriptive tools rather than portability cures.

## Headline results

Under **LD-block** cross-validation (~109k test variants):

| Analysis | Result |
|----------|--------|
| AF + LD + selection features | AUROC **0.627** for high-\(I^2\) failure |
| AF alone → AF+LD | **0.619 → 0.629** |
| vs coarse \(F_{ST}\) | Higher (paired ΔAUROC CI excludes 0) |
| Feature family permutation | Largest AUROC drop from **AF**, then **LD** |
| Sign-discordance endpoint | AUROC **0.77** |
| Score edits (Catalog PGS) | Frequency-aware filters (\(F_{ST}\)/MAF) reduce ancestry mean score separation; matched-mass Monte Carlo supports several modes |
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
