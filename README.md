# Biological Transferability Atlas

**Open-data prediction of cross-ancestry PRS/GWAS portability failure — and honest score-edit interventions — without requiring clinical phenotypes.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/python-3.11-green.svg)](bio-transfer-atlas/environment.yml)

This repository implements the **Biological Transferability Atlas (BTA)** / FAIRGEN-Open empirical stack: variant-level models that predict where polygenic scores and GWAS effects fail across ancestries, using **cross-ancestry effect concordance** (\(I^2\), sign discordance) as labels, then stress-testing Catalog PGS edits on **1000 Genomes GRCh38**.

> **What this is:** an open audit layer for transferability risk (AF/LD-driven structure, leakage-aware CV, matched intervention controls, Duffy/PAGE honesty checks).  
> **What this is not:** clinical risk prediction, fairness certification, or a claim that editing weights “fixes” phenotype \(R^2\) in every ancestry.

---

## Why it exists

Cross-ancestry PRS loss is often blamed on “ancestry” as a monolith. Recent work (e.g. Hu et al.; Harpak et al.) emphasizes **allele-frequency and LD structure**, trait-specific biology, and the limits of phenotype-based fairness metrics. BTA operationalizes that view on **fully open** data:

1. **Predict** variant-level portability failure from AF / LD / selection features under **LD-block** cross-validation.  
2. **Benchmark** against \(F_{ST}\) and population-distance encodings (not a fake peer contest vs trait-level \(r_g\)).  
3. **Intervene** on Catalog scores and measure ancestry **mean separation (MAD)** with matched-n / matched-mass Monte Carlo nulls.  
4. **Anchor** biology (Duffy/WBC) and external naming (PAGE vs internal Pan-UKB sensitivity).

Primary claim numbers use **LD-block** AUROC (AF_LD_SEL ≈ **0.627–0.629**), not leaky variant holdouts.

---

## Headline results (freeze)

| Result | Evidence |
|--------|----------|
| AF+LD predict high-\(I^2\) failure | LD-block AUROC ≈ **0.627**; AF perm. drop ≫ LD ≫ SEL |
| Beats coarse \(F_{ST}\); soft vs distance | Paired Δ vs FST CI excludes 0; vs POP_DISTANCE does **not** |
| Interventions = mean separation | fst/maf shrink MAD; risk filters often worsen; matched-mass \(p\) often low |
| Duffy positive control | AFR null-hom vs &lt;2 ΔWBC PGS ≈ **0.119** [0.024, 0.206] |
| PAGE external | **n_variants = 4384** after allele QC; β corr ≈ 0.01 |
| Trees ≥ GAT | GAT ≈ 0.619 vs XGB ≈ 0.623; attention is **descriptive** |

Machine-readable tables: [`bio-transfer-atlas/results/tables/`](bio-transfer-atlas/results/tables/) and the freeze pack [`bio-transfer-atlas/BTA_robustness_freeze_results_20260720/`](bio-transfer-atlas/BTA_robustness_freeze_results_20260720/).  
Narrative: [`bio-transfer-atlas/paper/`](bio-transfer-atlas/paper/) · [`bio-transfer-atlas/results/BTA_M4_M6_RESULTS_REPORT.md`](bio-transfer-atlas/results/BTA_M4_M6_RESULTS_REPORT.md).  
Methods vision (aspirational architecture): [`agents/FAIRGEN_Open_Methods_v1.md`](agents/FAIRGEN_Open_Methods_v1.md).

---

## Repository layout

```
fairness/                          # git root
├── agents/                        # FAIRGEN methods + research notes
├── bio-transfer-atlas/            # runnable atlas codebase
│   ├── configs/                   # experiment YAML
│   ├── scripts/                   # pipeline + evaluation entrypoints
│   ├── src/bta/                   # library code
│   ├── workflow/                  # Snakemake
│   ├── paper/                     # manuscript freeze
│   ├── results/tables/            # small CSV artifacts (tracked)
│   ├── tests/
│   ├── environment.yml
│   └── Makefile
└── LICENSE
```

**Not in git** (by design): genotypes (`.pgen`/VCF), large parquet lakes, model weights, raw Pan-UKB dumps. See `bio-transfer-atlas/data/**/README.md` / manifests for provenance. Heavy compute typically runs on a GPU host (e.g. Lambda) with project root mirrored under NFS.

---

## Quickstart

```bash
cd bio-transfer-atlas
mamba env create -f environment.yml
mamba activate bta

# Optional: smoke Snakemake on chr22 configs
make smoke

# Unit tests
make test

# Robustness gate (requires result tables present)
python scripts/gate_literature_roadmap.py
```

### Common evaluation scripts

| Script | Purpose |
|--------|---------|
| `scripts/run_nested_ablation_and_paired_auroc.py` | AF/LD/SEL nest + paired ΔAUROC verdicts |
| `scripts/run_grouped_permutation_importance.py` | Family permutation attribution |
| `scripts/evaluate_intervention.py` | MAD + retention; hooks matched MC |
| `scripts/run_intervention_matched_random_controls.py` | 500-draw matched-n/mass + LOSO |
| `scripts/run_duffy_positive_control.py` | Allele-audited Duffy / ACKR1 split |
| `scripts/run_external_page_grch38.py` | PAGE liftover + allele QC ladder |
| `scripts/gate_literature_roadmap.py` | Freeze gate for paper claims |

Download / preprocess scripts assume network access and substantial disk (see `LOCAL_COMPUTE.md`, `ROADMAP_TO_COMPLETION.md`).

---

## Claim discipline (read before citing)

1. Cite **LD-block** metrics as primary — not variant-holdout AUROC.  
2. Do **not** phrase results as “AF/LD beats Popcorn \(r_g\)”; trait-constant Z-concordance is a diagnostic, not a peer.  
3. Soften AF_LD_SEL vs **POP_DISTANCE** unless the paired CI excludes zero.  
4. Interventions reduce ancestry **mean separation**, not proven phenotype accuracy.  
5. SuSiE `fine_mapped` requires **signed** LD (`--r-unphased`); GAT attention ≠ mechanism.  
6. PAGE is external only with `n_variants ≥ 500` + QC counts; Pan-UKB |Δβ| filters are **`internal_only`**.

Run `python bio-transfer-atlas/scripts/gate_literature_roadmap.py` before treating the draft as frozen.

---

## Relation to FAIRGEN-Open

[`agents/FAIRGEN_Open_Methods_v1.md`](agents/FAIRGEN_Open_Methods_v1.md) describes the full theoretical stack (three-way AF/LD/selection decomposition, PopSpec/HPRN, pathway atlas, local-ancestry RFMix2). **This repo’s delivered paper path** is the open atlas + intervention bake-off: phenotype-free concordance labels, nested AF/LD/SEL evidence (AF-dominant), matched intervention nulls, and honesty controls. Heavier architecture pieces (PopSpec pretrain, full PRSM-PORT, RFMix2) remain roadmap unless marked complete in results reports.

---

## Data sources

| Source | Role |
|--------|------|
| 1000 Genomes / IGSR (GRCh38) | Genotypes for scoring & AF/LD features |
| PGS Catalog | Published weights |
| Pan-UK Biobank | Multi-ancestry sumstats → \(I^2\) / sign labels |
| PAGE (GWAS Catalog) | External β concordance after liftover |
| Reactome / constraint / selection tracks | Optional pathway & SEL features |

Respect each source’s license and terms; this repo’s **code** is MIT.

---

## Citation

```
Biological Transferability Atlas of Polygenic Scores Across Global Populations
[preprint DOI TBD]
```

If you use the FAIRGEN conceptual framing, also cite the methods note in `agents/FAIRGEN_Open_Methods_v1.md` and the primary biology papers it builds on (Hu et al.; Harpak et al.).

---

## License

Code: [MIT](LICENSE).  
Data: retained by original providers — see manifests under `bio-transfer-atlas/data/`.
