# Why polygenic scores fail across ancestry — and what actually helps

**Biological Transferability Atlas · M1–M3 campaign (July 2026)**  
Compute: Lambda A100 40GB · project root `/lambda/nfs/geeg/fairness`

> **Paper draft (M4–M6 freeze):** [`../paper/`](../paper/) (Methods / Results / Discussion).  
> **M4–M6 campaign write-up:** [`BTA_M4_M6_RESULTS_REPORT.md`](BTA_M4_M6_RESULTS_REPORT.md).  
> This file remains the M1–M3 ledger; manuscript claim language lives under `paper/`.

---

## One-sentence takeaway

Frequency and linkage differences between populations are better early warning signs of score failure than “how genetically distant” two groups are — and when we try to fix score shifts, ancestry-aware weight surgery beats simply deleting variants the model thinks are risky.

---

## 1. The intuition (start here)

A **polygenic score** is a weighted checklist of DNA variants used to estimate genetic risk or a trait. Most scores are built mainly from European-ancestry studies. When you apply the same checklist to other groups, two things often go wrong:

1. **The checklist misses the real drivers.** The variants that tag a disease-related region in one population may be different — or differently correlated — in another. That correlation structure is called **linkage disequilibrium (LD)**.
2. **The weights assume the wrong frequencies.** A variant that is rare in Europe can be common elsewhere (**allele frequency, AF**). That changes how much the score moves.

This project asks a practical question: can we **predict which variants will break transfer** before we rely on the score, and can we **edit the score** so ancestry-related shifts shrink without destroying the signal?

---

## 2. What we did in this campaign

Work ran on an A100 GPU machine with large RAM. We used the GPU for gradient-boosted trees and neural nets. Genotype scoring with PLINK stayed on CPU (that software does not accelerate on GPU), so we parallelized those jobs across chromosomes instead.

### Predict failure
Label variants that show high inconsistency across ancestries (heterogeneity / sign discordance), then train classifiers on different feature “stories”:

- AF + LD selection patterns  
- Simple genetic distance between populations  
- Proxy correlation features built from AF/LD summaries  
- A chromosome-22 estimate of **real genetic correlation** between ancestries  

We also tested whether **AlphaMissense** (protein-damage scores for missense variants) adds predictive value on top of AF/LD.

### Measure and intervene
Scored nine published polygenic scores on ~3,200 people from the 1000 Genomes panel across all 22 autosomes, then measured how far non-European group means sit from European means. Applied ten edit strategies (drop “risky” variants, reweight, FST/MAF/LD filters, random controls) and re-scored to see which edits shrink those shifts.

---

## 3. How far did we expand? (full autosome vs pilots)

“Genome-wide” here means the **22 autosomes** (chr1–22), not sex chromosomes. Coverage is uneven by design: the score atlas is autosome-complete; some new feature families are still chromosome-22 pilots.

| Analysis layer | What it covers | Status |
|----------------|----------------|--------|
| PGS scoring on 1000 Genomes | chr1–22 | Full autosome |
| Variant failure-prediction table | ~37.6M rows, chr1–22 | Full autosome |
| Score-shift atlas + interventions | 9 scores × 10 modes | Full autosome |
| Real genetic-correlation features | chr22 Z-correlations | Pilot |
| Fine-mapping credible sets (SuSiE) | chr22, 4 traits × 2 ancestries | Pilot |
| Within-ancestry AF heterogeneity | chr8–22 (not yet 1–7) | Partial |

---

## 4. Results that matter

### A. Predicting which variants will fail to transfer

We report **AUROC**: the chance that a random “fails to transfer” variant ranks above a random “transfers fine” variant. 0.5 is coin-flip; higher is better. The fair comparison uses **LD-block cross-validation** (train and test on different genomic neighborhoods) so the model cannot cheat by memorizing nearby correlated variants.

| Feature story | AUROC (LD-block) |
|---------------|------------------|
| **AF + LD selection** | **0.627** |
| RG proxy (AF/LD mix) | 0.618 |
| Population distance | 0.617 |
| Single FST | 0.580 |
| Real RG (chr22 pilot) | 0.566 |

**In plain language:** knowing how allele frequencies and LD patterns differ between populations is more useful, right now, than knowing a single “genetic distance” number or a chromosome-22 genetic-correlation estimate.

**AlphaMissense:** AF+LD alone ≈ AF+LD+AlphaMissense (~0.73 on the VEP-style eval). No meaningful lift — likely because few scored variants are missense with usable AlphaMissense coverage.

GPU models agreed: gradient-boosted trees reached ~0.62 AUROC on the same LD-block task; a large neural net reached ~0.61. Attribution pointed at **frequency features** as the dominant drivers. A looser variant-level train/test split produced much higher AUROC (~0.83) — that number is **not** the one to cite for claims, because nearby variants leak information across the split.

### B. Genome-wide scores and the portability model

| Quantity | Value |
|----------|-------|
| Samples × polygenic scores | 3,202 × 9 |
| Portability model AUROC (variant split) | 0.687 |
| Master variant×trait rows | ~37.6M |

Scores were summed over chr1–22 on the GRCh38 1000 Genomes panel. The portability model predicts high cross-ancestry inconsistency for individual variants.

### C. Do interventions shrink ancestry score shifts?

For each score we measure mean absolute distance of non-European ancestry groups from the European mean (after scaling by European spread). **Reduction** = baseline distance minus post-edit distance. Positive means the edit made groups more similar to Europe on that score; negative means the edit made disparities worse.

| Edit strategy | Mean MAD reduction (9 scores) |
|---------------|-------------------------------|
| **FST filter** | **+0.50** |
| MAF filter | +0.21 |
| Flag only | +0.05 |
| Reweight (exp / linear) | ~0 to +0.03 |
| Random drop (control) | +0.01 |
| LD filter / drop top-risk 5–20% | **negative** (often worse) |

The strongest single wins were FST-based edits on several coronary and type-2-diabetes scores (reductions around +0.7 to +1.2 on individual scores). Dropping the variants the portability model flagged as highest-risk often **increased** mean ancestry shifts — a caution against treating “predicted failure” as an automatic deletion list without checking score-level consequences.

---

## 5. What claims this supports (and what it does not)

| Claim | Verdict | Why |
|-------|---------|-----|
| AF/LD patterns predict transfer failure better than simple distance or chr22 RG | **Supported** (LD-block) | 0.627 vs 0.566–0.618; same ordering on GPU trees |
| AlphaMissense / coding constraint improves prediction here | **Not supported** | No AUROC lift once AF/LD features are present |
| Mechanisms are largely frequency-driven | **Supported** | Attribution dominated by AF features |
| We can reduce ancestry score shifts with targeted edits | **Supported, mode-dependent** | FST/MAF help on average; risk filters often worsen MAD |
| Genome-wide genetic correlation is useless | **Do not claim** | RG_REAL is still a chr22 pilot, not autosome Popcorn \(r_g\) |
| Fine-mapping fully explains failure genome-wide | **Do not claim** | SuSiE outputs are a chr22 pilot with fallbacks |

---

## 6. Why it matters

Equity in genomic medicine is not only about collecting more diverse GWAS — it is also about knowing **when an existing score is lying** and **how to repair it** without waiting years for a perfect multi-ancestry rebuild.

- **For method developers:** prioritize frequency and LD diagnostics over coarse ancestry-distance summaries when screening variants for portability risk.
- **For score curators:** ancestry-aware filters (differentiation / frequency) look more promising than blunt “delete the model’s scary variants” policies — at least for reducing mean score shifts on this panel.
- **For the next milestone:** finish subpop AF on chr1–7, and only then decide whether to push real genetic correlation and fine-mapping from chr22 pilots to full autosome — those are the pieces that still limit how strongly we can word mechanism claims.

---

## Primary artifacts

- `results/LAMBDA_M1_M3_NOTES.md`
- `results/tables/ablation_ldblock_and_baselines_genomewide.csv`
- `results/tables/ablation_xgboost_gpu_expanded.csv`
- `results/tables/intervention_results.genomewide.csv`
- `results/tables/intervention_summary.genomewide.txt`
