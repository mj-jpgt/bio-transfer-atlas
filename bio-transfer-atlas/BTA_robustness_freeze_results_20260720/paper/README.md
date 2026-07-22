# Biological Transferability Atlas

Manuscript draft freeze (M6). Numbers reflect genome-wide LD-block analyses unless noted.
Update AUROCs after autosomal \(r_g\) re-ablation completes on Lambda.

## Abstract (draft)

Cross-ancestry polygenic scores often lose accuracy outside European discovery cohorts.
We build an open-data **Biological Transferability Atlas**: variant-level models that predict
GWAS/PRS portability failure from allele-frequency (AF), LD, and selection features, then
test score-edit interventions on 1000 Genomes GRCh38 genotypes. Under LD-block
cross-validation, AF+LD features achieve AUROC ≈ **0.627**, exceeding coarse \(F_{ST}\) and
chr22 genetic-correlation proxies; population-distance features are competitive (~0.617)
and should not be dismissed as uniformly weak. Frequency-aware interventions (\(F_{ST}\)/MAF
filters) reduce mean absolute EUR–non-EUR score gaps on average, whereas pruning top
predicted-risk variants often worsens gaps. Positive controls (Duffy/WBC; autoimmune MHC
stress-tests) and open external sumstat concordance with bootstrap CIs anchor biological
claim discipline. Fine-mapping tiers and an LD-graph GAT provide descriptive attention weights without
claiming that SuSiE or graph nets “solve” allele-frequency bias.

See [INTRODUCTION.md](INTRODUCTION.md), [METHODS.md](METHODS.md), [RESULTS.md](RESULTS.md), [DISCUSSION.md](DISCUSSION.md).
