# Introduction

Polygenic scores (PGS) and GWAS effect estimates routinely lose accuracy outside the
ancestries in which they were discovered. Allele-frequency (AF) differences and linkage-
disequilibrium (LD) tagging changes are leading mechanistic explanations; population
genetic distance and trait-level genetic correlation (\(r_g\)) are often cited as alternatives
or complements. Fine-mapping and graph neural networks have been proposed as ways to
isolate “causal” signal, but whether they overturn AF-driven failure genome-wide remains
an open empirical question.

We construct an open-data **Biological Transferability Atlas** that (i) predicts **variant-
level** portability failure from AF/LD/selection features under **LD-block** cross-validation,
(ii) benchmarks those predictors against \(F_{ST}\) and richer population-distance encodings
(trait-level Z-score concordance is evaluated at the **trait scale**, not as a peer of AF/LD
in variant ablation), and (iii) stress-tests score-edit interventions—including a Duffy /
WBC positive control and MHC-heavy autoimmune scores—on 1000 Genomes GRCh38 genotypes.
Interventions are scored on ancestry **mean separation**, not phenotype accuracy.
Open sumstat analyses separate **internal Pan-UKB concordance sensitivity** from
**external PAGE** validation after GRCh38 liftover (n≥500 required).

Our stance follows recent synthesis (Wang/Dahl and related work): AF and LD remain
primary; fine-mapping tiers (signed-LD SuSiE only) and GAT attention are explanatory /
descriptive lenses, not claimed cures for portability gaps. Claim numbers in Results use
LD-block AUROC (≈0.627 for AF_LD_SEL) rather than leaky variant holdouts.
