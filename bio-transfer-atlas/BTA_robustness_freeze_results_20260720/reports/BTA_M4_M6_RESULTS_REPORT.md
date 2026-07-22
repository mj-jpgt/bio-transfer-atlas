# Paper freeze write-up — robustness fixes (M4–M6 revised)

**Biological Transferability Atlas · robustness pass before freeze (July 2026)**  
Compute: Lambda A100 · project root `/lambda/nfs/geeg/fairness`

Manuscript: [`../paper/`](../paper/)  
Earlier ledger: [`BTA_M1_M3_RESULTS_REPORT.md`](BTA_M1_M3_RESULTS_REPORT.md)

---

## One-sentence takeaway

Under LD-block CV, AF+LD features predict variant-level high-\(I^2\) failure; trait-level Z-score concordance is **not** a fair peer at that scale; interventions reduce ancestry **mean separation** (not proven accuracy); SuSiE primary claims require signed LD; PAGE external validation requires GRCh38 join with n≥500.

---

## Variant-scale peer contest (LD-block CV)

| Feature group | AUROC |
|---------------|------:|
| **AF_LD_SEL** | **0.629** |
| AF_LD | 0.627 |
| AF | 0.616 |
| POP_DISTANCE | 0.619 |
| FST | 0.605 |

`TRAIT_CONSTANT_Z_DIAGNOSTIC` (formerly RG_REAL) is **excluded** from the peer contest and from AF nest features (no FST-in-AF). Phrase: *trait-level concordance alone provided little discrimination among variant-level heterogeneity labels*—not “AF/LD beats \(r_g\)”.

Paired ΔAUROC (incl. AF_LD_SEL−AF, AF_LD−AF, AF_LD_SEL−AF_LD) with verdicts: `auroc_paired_delta_ldblock.csv`. Nested AF/LD/SEL: `ablation_nested_af_ld_sel.csv`. Grouped permutation: `grouped_permutation_importance_ldblock.csv` (AF drop ≫ LD ≫ SEL). Soften vs POP_DISTANCE (`no_clear_difference`).

Pearson Z correlations are labeled **`panukbb_z_concordance_*`** / estimand `cross_ancestry_z_score_concordance`—not Popcorn genetic-effect \(r_g\).

### Trait-scale table

`trait_scale_portability.csv`: one row per trait×ancestry-pair with Z-concordance, n variants, trait class, mean high-\(I^2\) rate. Concordance is a legitimate **aggregate** predictor here.

---

## Concordance analyses (naming)

### Internal Pan-UKB sensitivity (not external)

`internal_panukbb_concordance_sensitivity.csv` — LDL EUR–AFR chr22:

| Metric | Value |
|--------|------:|
| β correlation (raw) | −0.005 [−0.011, −0.000] |
| After drop top 10% \|Δβ\| | 0.118 [0.110, 0.126] |
| Matched random drop | ≈ unchanged (null) |

Improvement after dropping largest \|Δβ\| is **expected by construction**. Status: `internal_only`.

### External PAGE (GRCh38)

`external_page_validation.csv` — liftover hg19→hg38, allele QC ladder, filter frozen on Pan-UKB only, LD-block bootstrap. See also `external_page_qc_counts.csv`.

| Metric | Value |
|--------|------:|
| n_variants after allele QC | **4384** (`ok_external`) |
| input → final (qc) | 5416 → 4384 (800 ambiguous, 181 indel, 51 unresolved) |
| β correlation (Pearson) | 0.010 [−0.027, 0.042] |
| Spearman / sign concordance | −0.022 / 0.49 |
| After Pan-UKB-frozen filter | 0.013; matched random keep does not help |

Held-out-trait ranking companion: `external_heldout_trait_ranking.csv`.

---

## SuSiE / fine-mapping

Primary `fine_mapped` requires `fallback_mode == signed_ld` from `plink2 --r-unphased square`
(PLINK2 a.7 writes `.unphased.vcor1`). PSD ridge; `estimate_residual_variance=FALSE`.
Unsigned \|r\| / identity / \|z\| / legacy `real_ld` → `pipeline_fallback`, non-primary.
Meta: `data/labels/susie/susie_primary_meta.json`. CAD EUR/AFR currently majority
`signed_ld`; BMI/LDL/T2D continue on the signed-LD / vcor1 path.

---

## Interventions = ancestry mean separation

Primary metric remains MAD (mean \|Δ EUR\| score gap)—**not** phenotype portability/accuracy. Retention / variance / corr(edited, original): `intervention_retention_variance_metrics.csv`. Matched-n / matched-mass Monte Carlo (≥500 draws) + empirical \(p\): `intervention_matched_random_controls.csv`. LOSO: `intervention_loso_mad_by_mode.csv`. Curve: `fig_intervention_retention_curve.png`.

| Clade | Best mode (MAD) | Notes |
|-------|-----------------|-------|
| Metabolic | fst / maf | Risk filters often worse |
| WBC | maf | Mean separation |
| Autoimmune | fst | MHC-heavy |

---

## Duffy positive control (allele-audited)

rs2814778 GRCh38 chr1:159204893 **T/C**; PLINK counted **T**; dosage flipped so **`duffy_null_dose` = copies of C** (Duffy-null). AFR null AF ≈ 0.96; EUR ≈ 0.006 (`duffy_allele_audit.csv`, `duffy_null_af_by_superpop.csv`).

Primary AFR contrast: **dose==2 vs dose&lt;2** with bootstrap CIs. ACKR1 score decomposition (lead SNP / window / outside) → `duffy_ackr1_score_decomposition.csv` (`ackr1_fraction`).

| Ancestry | Null copies | n | Mean WBC PGS | Underpowered |
|----------|-------------|--:|-------------:|:------------:|
| AFR | 0 | 4 | 14.95 | yes |
| AFR | 1 | 40 | 15.15 | no |
| AFR | 2 | 617 | 15.26 | no |
| EUR | 0 | 497 | 14.74 | no |
| EUR | ≥1 | 6 | 14.62 | yes |

Do **not** claim a linear AFR dose-response from the n=4 cell; prefer homozygous / dominant contrasts. Caption must list alleles and counts.

---

## GAT / mechanism language

GAT AUROC ≈ 0.619 vs XGBoost ≈ 0.623. Trees primary. Attention = **descriptive graph weighting** / `descriptive_attention_weights`, not mechanism evidence. Prefer SHAP / grouped permutation (`grouped_permutation_importance_ldblock.csv`).

---

## Interpretation notes

Primary mechanism claim is LD-block AF+LD+SEL (~0.629 AUROC), not variant-holdout.
Interventions report ancestry mean separation (MAD), not phenotype accuracy.
SuSiE primary uses signed LD only; Duffy uses homozygous 2-vs-&lt;2 plus ACKR1 fraction;
PAGE external requires allele QC counts with adequate variant overlap.
