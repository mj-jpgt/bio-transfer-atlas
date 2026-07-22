# Discussion

## What transfers—and what does not

Allele-frequency and LD structure remain the dominant, portable signals of GWAS/PRS
failure across ancestries, aligning with Wang/Dahl-style accounts of AF-driven
heterogeneity. Our contribution is operational: a **variant-level risk model** plus an
**intervention bake-off** on open genotypes, with claim discipline under LD-block CV.

Coarse \(F_{ST}\) underperforms AF+LD (paired ΔAUROC CI excludes zero), but
**population-distance** features nearly match AF_LD_SEL (~0.619 vs ~0.627); the paired
ΔAUROC CI includes zero, so we do not claim a large AF/LD vs distance gap. Nested AF /
LD / SEL ablations show AF+LD carries most of the signal.

Trait-level Z-score concordance is a **trait-scale** predictor, not a fair peer of
variant-level AF/LD. We therefore do not phrase results as an AF/LD vs \(r_g\) peer contest.

## Interventions reduce mean separation, not proven accuracy

Frequency-aware edits (\(F_{ST}\)/MAF) shrink mean absolute ancestry gaps on average.
Aggressively dropping high predicted-risk variants often **increases** gaps. Retention,
variance, and correlation metrics—and matched random controls—separate compression from
biology. We do **not** claim improved phenotypic portability without phenotype association.

Duffy gating is a biology-informed special case, not a general recipe.

## Fine-mapping and graphs as lenses, not cures

SuSiE credible sets (signed LD only) and GAT attention clarify *where* portability risk
concentrates (e.g., ACKR1 for WBC) but do not overturn AF dominance genome-wide. We refuse
silent \(R=I\) or unsigned-|r| “SuSiE” claims. GAT attention is descriptive graph
weighting—not mechanism evidence.

## Concordance scope

Internal Pan-UKB |Δβ| filters are sensitivity checks (expected-by-construction). External
PAGE after GRCh38 liftover joins n=5416 variants with a filter frozen without looking at
PAGE; concordance remains near zero and the frozen filter does not beat a matched random
keep-set. Phenotype-level multi-ancestry \(R^2\) remains future work.

## Limitations

- Labels derive from Pan-UKB stratified GWAS; Catalog PGS discovery ancestries vary.
- Z-concordance is not Popcorn genetic-effect \(r_g\) until `popcorn compute`+`fit` succeeds.
- SuSiE signed-LD coverage depends on PLINK `--r-unphased` extract success per block.
- 1000G is a reference panel, not a phenotype cohort; MAD is mean separation, not clinical
  calibration.
- Autoimmune MHC structure may dominate RA/IBD behavior beyond AF features.
- Duffy dose-2 cells can be tiny (e.g. AFR n≈4); avoid overclaiming linear dose-response.

## Claim checklist (freeze)

1. Cite **LD-block** AF_LD_SEL AUROC (~0.627), not leaky variant holdout, as primary.
2. Do **not** pit trait-constant Z-concordance against AF/LD as a peer contest.
3. Soften AF vs POP_DISTANCE unless paired ΔAUROC CI excludes 0.
4. State interventions reduce ancestry **mean separation**; risk pruning can backfire.
5. AlphaMissense null + sparse-overlap caveat.
6. SuSiE primary = signed LD only; GAT descriptive vs trees.
7. Duffy allele audit + CIs; PAGE external only if n≥500 else `internal_only`.

Point machine-readable M1–M3 notes at this folder:
[`../results/BTA_M1_M3_RESULTS_REPORT.md`](../results/BTA_M1_M3_RESULTS_REPORT.md).
