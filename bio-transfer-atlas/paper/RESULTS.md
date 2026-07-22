# Results

## Variant-scale portability prediction (LD-block CV)

Under LD-block cross-validation, **AF_LD_SEL** predicts high-\(I^2\) failure with AUROC ≈
**0.629**. Coarse **FST** is weaker (≈0.605; paired ΔAUROC CI excludes zero; verdict
`higher`). **POP_DISTANCE** reaches ≈ **0.619**—close to AF+LD; paired ΔAUROC ≈ 0.010
[−0.0005, 0.018] does **not** exclude zero (`no_clear_difference`). Nested AF / LD / SEL
(`ablation_nested_af_ld_sel.csv`; AF excludes FST): AF alone ≈0.616; AF+LD ≈0.627;
AF_LD_SEL ≈0.629; SEL adds little beyond AF/LD (AF_LD_SEL−AF_LD CI includes 0). Grouped
permutation importance (`grouped_permutation_importance_ldblock.csv`): largest AUROC drop
from permuting AF, then LD; SEL/distance/FST near zero.

Variant-holdout AUROCs (~0.73) are **not** the primary claim (LD leakage).

**Trait-constant Z-concordance is not a peer comparator at the variant scale.** Joining a
four-value (trait×pair) concordance feature into a variant classifier only exploits
between-trait failure-rate differences. We report it only as
`TRAIT_CONSTANT_Z_DIAGNOSTIC` and phrase the result as: *trait-level Z-concordance alone
provided little discrimination among variant-level heterogeneity labels*—not an AF/LD vs
\(r_g\) peer contest.

Pearson correlation of ancestry Z-scores is labeled **cross-ancestry Z-score concordance**
(`panukbb_z_concordance_*`), not Popcorn genetic-effect \(r_g\).

## Trait / population-pair scale

At the **trait×ancestry-pair** scale, Z-concordance (and eventually true Popcorn \(r_g\)),
distance summaries, trait class, and GWAS size are legitimate predictors of aggregate
outcomes (`trait_scale_portability.csv`).

## Fine-mapping tiers

Primary `fine_mapped` requires SuSiE credible sets from **signed** LD
(`plink2 --r-unphased square`). Identity-\(R\), |z|-weights, and unsigned |r| from r² are
`pipeline_fallback` / non-primary. Stance: explanatory lens; AF differences still expected
to dominate.

## LD-graph GAT

GAT AUROC ≈ 0.619 vs XGBoost ≈ 0.623. Trees remain primary. Attention weights are
**descriptive graph visualizations**, not mechanism evidence.

## Interventions (ancestry mean separation)

Primary outcome is **mean absolute EUR–non-EUR score separation (MAD)**—not phenotype
accuracy or “repaired portability.” Frequency-aware edits (fst/maf) reduce mean separation
more than risk filters on average; report retention of variants, |β| mass, within-ancestry
score variance, and corr(edited, original) (`intervention_retention_variance_metrics.csv`).
Matched-n / matched-mass Monte Carlo (≥500 draws; AF-expected MAD) yields
`delta_vs_random_n`, `delta_vs_random_mass`, and empirical \(p\) vs random compression
(`intervention_matched_random_controls.csv`). Frequency-aware edits often beat
matched-mass random (low `empirical_p_mass`) but not matched-*n* random (high
`empirical_p_n` when mass differs). Leave-one-score-out: `intervention_loso_mad_by_mode.csv`.
Retention–MAD curve: `fig_intervention_retention_curve.png`.
| Clade | Best mode (MAD) | Notes |
|-------|-----------------|-------|
| Metabolic | fst / maf | Risk filters often worse |
| WBC | maf ≈ 0.345 | Mean separation, not accuracy |
| Autoimmune | fst ≈ 0.439 | MHC-heavy |

## Duffy positive control (WBC)

rs2814778 allele audit (`duffy_allele_audit.csv`): REF/ALT, PLINK counted allele, flip to
**Duffy-null allele C**, AF by superpop. Primary AFR biology contrast: **null-homozygous (dose==2) vs dose&lt;2** with bootstrap
CIs—AFR Δmean WBC PGS ≈ **0.119** [0.024, 0.206] (n_hom2=617, n_lt2=44). Keep 0/1/2
strata; n&lt;10 underpowered. ACKR1 score decomposition (PGS000191): window (3 weights)
Δ ≈ +0.57 vs outside Δ ≈ −0.45; `ackr1_fraction = Δ_window/Δ_total` can exceed 1 when
components oppose (`duffy_ackr1_score_decomposition.csv`). Lead SNP rs2814778 itself is
absent from this Catalog score (n_snp_weights=0).

## Concordance analyses (naming)

**Internal Pan-UKB sensitivity** (`internal_panukbb_concordance_sensitivity.csv`): EUR–AFR
β concordance before/after dropping top |Δβ|; improvement is expected by construction;
compare to matched random drop. **Not external validation.**

**External PAGE** (`external_page_validation.csv`): GRCh38 liftover, allele QC ladder
(match/swap/complement; drop ambiguous A/T·C/G and indels), filter frozen on Pan-UKB only
(`assert` PAGE columns disjoint from filter training). After QC: **n_variants=4384**
(`ok_external`; 800 ambiguous + 181 indels + 51 unresolved removed from 5416 position
joins; `external_page_qc_counts.csv`). Pearson ≈ 0.010; Spearman ≈ −0.022; sign
concordance ≈ 0.49. Held-out-trait ranking companion: `external_heldout_trait_ranking.csv`.

## Subpopulation AF and robustness

Subpop AF features cover chr1–22. Nested ablations and paired ΔAUROC CIs support claim
discipline around AF/LD vs distance.
