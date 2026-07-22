# Roadmap to Completion — Biological Transferability Atlas

A task-by-task checklist to take the project from the current **chr22 end-to-end
proof** to a **genome-wide, multi-source, publication-ready** atlas. Every task
includes (a) the **resources/links** you need, (b) a **copy-paste agent prompt**,
and (c) a **validation harness** (the acceptance gate that proves the task is done).

> Convention: prompts assume system Python + bundled `plink2`, `curl.exe`, and the
> existing `bio-transfer-atlas/scripts/` layout. Run gates with
> `python -m pytest tests/ -k <gate>` or the standalone gate scripts in `tests/`.

---

## 0. Status snapshot (already complete)

- [x] **Milestone A — Score-shift atlas** (`results/tables/score_shifts_*.parquet`,
      `results/figures/fig_atlas_heatmap.png`, `fig_distance_sensitivity.png`,
      `fig_rank_instability.png`)
- [x] **Stage 8 — Pan-UKBB chr22 labels** (`data/labels/gwas_concordance_labels.parquet`,
      490k variant-trait rows, 5–6 ancestries; tabix-range fetcher in
      `scripts/download_panukbb_chr22.py`)
- [x] **Stage 9 — Selection-turnover proxies** (PBS + Hudson F_ST,
      `data/features/selection/selection_turnover_features.parquet`)
- [x] **Stage 10 — Master table + leakage-safe splits**
      (`data/modeling/master_variant_table.parquet`)
- [x] **Stage 11 — Mechanism-ablation baselines** (`results/tables/ablation_*.csv`):
      AF+LD AUROC **0.87** vs F_ST **0.61** (unseen-variant split); trait-holdout ≈ chance.
- [x] **Stage 12 — Pathway-risk atlas** (`results/tables/pathway_risk_table.parquet`)

**Known gaps to close:** single-biobank non-EUR power (median I²≈0), no cross-trait
generalization on chr22, chr22-only, proxy-only selection features, no neural model,
no formal CIs / negative controls / reproducible DAG.

---

## Master resource table (verified links)

| Resource | Use | Link | Build |
|---|---|---|---|
| **Pan-UKBB** flat files + tabix | multi-ancestry GWAS betas | manifest: `https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz` ; files under `https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/` | GRCh37 |
| **BioBank Japan (BBJ)** | independent **EAS** GWAS | downloads: https://pheweb.jp/downloads · NBDC hum0197: https://humandbs.dbcls.jp/en/hum0197 · file dict: https://humandbs.dbcls.jp/files/hum0197/Dictfile_BBJ.html | GRCh37 |
| **FinnGen DF13** | independent **EUR (Finnish)** GWAS | access form: https://elomake.helsinki.fi/lomakkeet/124935/lomake.html · docs: https://finngen.gitbook.io/documentation/data-download · R12 manifest (public GCS): `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv` | GRCh38 |
| **FinnGen + Pan-UKBB meta** | ready-made cross-biobank meta (867 phenos) | https://public-metaresults-fg-ukbb.finngen.fi/ | GRCh38 |
| **gnomAD v4.1 constraint** | LOEUF / pLI / missense-Z (selection) | `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv` | GRCh38 |
| **UCSC phyloP/phastCons** | conservation (selection) | hg38 bigWig: `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw` | GRCh38 |
| **Ensembl 110 GFF3** | gene coordinates (per-chrom) | `https://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens/Homo_sapiens.GRCh38.110.chromosome.<CHR>.gff3.gz` | GRCh38 |
| **Reactome** | pathway hierarchy | https://reactome.org/download-data (`Ensembl2Reactome_All_Levels.txt`, `ReactomePathwaysRelation.txt`) | n/a |
| **PGS Catalog** | scoring files | https://www.pgscatalog.org/ · API: https://www.pgscatalog.org/rest/ | hmPOS GRCh37/38 |
| **1000G GRCh38** | reference panel / AF / LD / PCA | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/` | GRCh38 |
| **liftOver chains** | build conversion | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz` (+ reverse) | n/a |
| **selscan** | iHS / nSL / XP-EHH (real selection) | https://github.com/szpiech/selscan | n/a |

---

## Phase 1 — Strengthen labels with independent sources (highest priority)

> Why first: this directly addresses the median-I²≈0 / trait-holdout-≈-chance
> weakness. Better-powered independent EAS (BBJ) and EUR (FinnGen) estimates turn
> noisy single-biobank heterogeneity into real cross-source concordance signal.

### 1.1 BioBank Japan (EAS) chr22

- [x] Download BBJ chr22 for T2D, CAD, BMI, LDL; harmonize to GRCh38 targets.

**Resources:** https://pheweb.jp/downloads (per-trait `wget` list) · column dict:
https://humandbs.dbcls.jp/files/hum0197/Dictfile_BBJ.html (GRCh37, columns
`SNP, CHR, POS, REF, ALT, Frq, BETA, SE, P`). BBJ trait names: *Type 2 diabetes*,
*Coronary artery disease* (or *Ischemic heart disease*), *Body mass index*,
*LDL cholesterol*.

**Agent prompt:**
```
Create scripts/download_bbj_chr22.py. BBJ flat files are GRCh37, gzip (not bgzf-tabix),
~hundreds of MB each, genome-wide. For T2D/CAD/BMI/LDL: stream-download each
GWASsummary_<TRAIT>_Japanese_SakaueKanai2020.auto.txt.gz from pheweb.jp, decompress
on the fly, keep only CHR==22 rows, and save data/raw/bbj/chr22/<trait>.chr22.parquet
with columns [chr,pos,ref,alt,beta,se,pval,af]. Delete the full download after filtering.
Reuse the curl_full helper pattern from download_panukbb_chr22.py. Print kept-row counts.
```

**Validation harness** (`tests/gate_bbj.py`):
```python
import pandas as pd, glob, sys
files = glob.glob("data/raw/bbj/chr22/*.chr22.parquet")
assert len(files) == 4, f"expected 4 BBJ traits, got {len(files)}"
for f in files:
    d = pd.read_parquet(f)
    assert (d["chr"].astype(str) == "22").all(), f"{f}: non-chr22 rows"
    assert d["beta"].notna().mean() > 0.9 and len(d) > 50_000, f"{f}: too few usable betas"
print("GATE PASS: BBJ chr22 downloaded & filtered")
```

### 1.2 FinnGen (EUR-Finnish) chr22

- [x] Register, download FinnGen DF13 chr22 for the 4 traits (already GRCh38 → no liftover). *(Implemented T2D/CAD endpoints directly from manifest; BMI/LDL kept as meta-only per FinnGen scope.)*

**Resources:** form https://elomake.helsinki.fi/lomakkeet/124935/lomake.html →
emailed GCS instructions. Manifest pattern (public bucket):
`https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv`
(swap `r12`→`r13` once DF13 manifest is published). FinnGen endpoints:
`T2D = E4_DM2`, `CAD = I9_CHD` / `I9_IHD`, `BMI` (use a meta source — FinnGen is
disease-focused, so for BMI/LDL prefer the **FinnGen+Pan-UKBB meta** browser).

**Agent prompt:**
```
Create scripts/download_finngen_chr22.py. Read the FinnGen R13 manifest TSV, select
endpoints E4_DM2 (T2D) and I9_CHD (CAD). FinnGen sumstats are bgzipped + tabix-indexed
and GRCh38. Reuse the tabix-range chr22 logic from download_panukbb_chr22.py
(parse .tbi linear index -> chr22 coffset -> HTTP range to EOF -> filter chr==22).
Save data/raw/finngen/chr22/<trait>.chr22.parquet with [chr,pos,ref,alt,beta,sebeta,pval,af_alt].
For BMI/LDL (not native FinnGen endpoints) skip and note them as meta-only.
Document any manual GCS auth needed in a top-of-file comment.
```

**Validation harness** (`tests/gate_finngen.py`):
```python
import pandas as pd, glob
files = glob.glob("data/raw/finngen/chr22/*.chr22.parquet")
assert len(files) >= 2, "expected >=2 FinnGen disease traits (T2D, CAD)"
for f in files:
    d = pd.read_parquet(f)
    assert (d["chr"].astype(str).isin(["22","chr22"])).all()
    assert len(d) > 50_000
print("GATE PASS: FinnGen chr22 ready (GRCh38, no liftover)")
```

### 1.3 Multi-source concordance labels

- [x] Recompute concordance using Pan-UKBB ancestries **+ BBJ EAS + FinnGen EUR** as
      independent sources; add a **source-holdout** evaluation axis.

**Agent prompt:**
```
Extend scripts/build_concordance_labels.py into build_concordance_labels_multisource.py.
Treat each (biobank, ancestry) as an independent effect estimate for a variant/trait:
Pan-UKBB {AFR,AMR,CSA,EAS,EUR,MID}, BBJ {EAS}, FinnGen {EUR-FIN}. Align all to GRCh38
target variant_ids (BBJ needs hg19->hg38 liftover via the cached map; FinnGen is already
hg38). Recompute sign_concordance, Cochran's Q, I2, het_pval, risk_class, plus a new
column n_sources and a per-source-pair sign-agreement matrix. Add source_group column so
downstream can do source-holdout (e.g., train on Pan-UKBB-only, test on BBJ/FinnGen).
Write data/labels/gwas_concordance_labels_multisource.parquet and refresh label_audit.txt.
```

**Validation harness** (`tests/gate_multisource_labels.py`):
```python
import pandas as pd
d = pd.read_parquet("data/labels/gwas_concordance_labels_multisource.parquet")
assert d["n_sources"].max() >= 3, "need >=3 independent sources for >=some variants"
assoc = d[d["associated"]]
# independent sources should reduce noise -> associated median I2 should be > 0
assert assoc["I2"].median() > 0.0, "multi-source did not improve heterogeneity signal"
print("GATE PASS: multi-source labels (median assoc I2 =", round(assoc['I2'].median(),3), ")")
```

---

## Phase 2 — Real mechanism features (replace proxies)

### 2.1 gnomAD constraint (LOEUF / pLI / missense-Z)

- [x] Add gene-level constraint to each variant (via gene mapping already built).

**Resource:** `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv`
(transcript-level; pick MANE Select else canonical; fields `lof.oe_ci.upper` = LOEUF,
`lof.pLI`, `mis.z_score`).

**Agent prompt:**
```
Create scripts/compute_constraint_features.py. Download gnomad.v4.1.constraint_metrics.tsv,
reduce to gene-level (MANE Select > canonical > longest CDS), keep ensg, LOEUF (lof.oe_ci.upper),
lof.pLI, mis.z_score. Join via data/annotations/variant_to_gene.parquet (nearest/overlapping
gene). For multi-gene variants take the most-constrained (min LOEUF). Save
data/features/selection/constraint_features.parquet keyed by variant_id. Add these columns
to the SEL feature group in feature_groups.json.
```

**Harness** (`tests/gate_constraint.py`):
```python
import pandas as pd
c = pd.read_parquet("data/features/selection/constraint_features.parquet")
assert {"variant_id","LOEUF","pLI","mis_z"} <= set(c.columns)
assert c["LOEUF"].notna().mean() > 0.5
print("GATE PASS: constraint features attached")
```

### 2.2 Conservation (phyloP) + recombination rate

- [ ] Add per-variant phyloP (and optionally recombination rate) features.

**Resources:** phyloP bigWig `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw`
(read with `pyBigWig`); recombination map (deCODE/1000G genetic maps).

**Agent prompt:**
```
Create scripts/compute_conservation_features.py. pip install pyBigWig. For each chr22 target
variant, query hg38.phyloP100way.bw at its position (and mean over +/-25bp). Save
data/features/selection/conservation_features.parquet [variant_id, phyloP, phyloP_win].
Add to SEL group. If pyBigWig is unavailable on Windows, document the WSL/conda fallback.
```

**Harness:** assert `phyloP` non-null > 90%, range roughly [-20, 10].

### 2.3 Real selection scans (iHS / nSL) — optional, strongest signal

- [ ] Run selscan on 1000G phased chr22 per superpop to get iHS/nSL.

**Resource:** https://github.com/szpiech/selscan ; 1000G phased VCFs (table above).

**Agent prompt:**
```
Create scripts/compute_ihs_features.py. Using bundled plink2 + selscan, for each superpop
(AFR,EUR,EAS,SAS,AMR): subset 1000G phased chr22 to that superpop, run `selscan --ihs`
and `--nsl` with the genetic map, normalize with norm, then attach standardized |iHS| and
nSL per target variant_id. Save data/features/selection/ihs_features.parquet. This replaces
PBS as the primary selection signal; keep PBS as a fallback.
```

**Harness:** assert per-superpop iHS columns exist and >70% variants have a value.

---

## Phase 3 — Genome-wide scaling

- [ ] Parameterize every chr22 script over `--chrom {1..22}` and run the full genome.

**Agent prompt:**
```
Refactor compute_af_features.py, compute_ld_features.py, compute_selection_features.py,
download_panukbb_chr22.py (-> download_panukbb_chrom.py), build_concordance_labels.py,
build_master_table.py, build_pathway_risk.py to accept --chrom and write per-chrom outputs,
then a concat step. For Pan-UKBB, generalize the tbi parser to any contig name. Keep memory
bounded by processing one chromosome at a time. Provide a driver scripts/run_genomewide.py
that loops 1..22 and concatenates into data/modeling/master_variant_table_genomewide.parquet.
Re-run Stage 11 ablation genome-wide and refresh results/tables + figures.
```

**Harness** (`tests/gate_genomewide.py`):
```python
import pandas as pd
m = pd.read_parquet("data/modeling/master_variant_table_genomewide.parquet")
assert m["variant_id"].str.split(":").str[0].nunique() >= 20, "need >=20 chromosomes"
assert len(m) > 1_000_000
print("GATE PASS: genome-wide master table", len(m), "rows")
```

---

## Phase 4 — HPRN-lite model (the architecture from the methods doc)

> Only build after Phase 1 gate passes (a mechanism group must beat F_ST under
> source/trait holdout). Otherwise baselines are the result.

- [ ] Implement a compact version of the Hierarchical Portability Risk Network:
      PopSpec encoder → LD-graph attention → structured mechanism bottleneck
      (AF / LD / selection heads) → Reactome-hierarchical pathway pooling.

**Resources:** `torch`, `torch_geometric` (GAT), Reactome relation file for pooling.
Reference: `FAIRGEN_Open_Methods_v1.md` §architecture.

**Agent prompt:**
```
Create scripts/train_hprn.py (PyTorch). Inputs: master_variant_table features grouped into
AF/LD/SEL blocks + LD adjacency (from compute_ld_features). Model: per-block encoders ->
a structured bottleneck with one latent per mechanism (AF, LD, selection) regularized to be
predictive of that block's variance -> Reactome-hierarchical pooling using
ReactomePathwaysRelation.txt -> heads for y_high_I2 (BCE), I2 (MSE), risk_class (CE).
Use split_variant and split_trait and source-holdout. Log AUROC/AUPRC/Spearman per head and
compare to the HGB baseline in ablation_classification.csv. Save results/tables/hprn_metrics.csv
and a mechanism-attribution table (which latent explains each prediction).
```

**Harness** (`tests/gate_hprn.py`):
```python
import pandas as pd
h = pd.read_parquet("results/tables/hprn_metrics.csv") if False else pd.read_csv("results/tables/hprn_metrics.csv")
base = pd.read_csv("results/tables/ablation_classification.csv")
best_base = base[(base.subset=="associated")&(base.split=="split_variant")&(base.model=="hgb")]["AUROC"].max()
assert h["AUROC"].max() >= best_base - 0.02, "HPRN must be within 2pts of / beat best baseline"
print("GATE PASS: HPRN trained; best AUROC", round(h['AUROC'].max(),3), "vs baseline", round(best_base,3))
```

---

## Phase 5 — Evaluation suite & statistical rigor

- [x] Bootstrap 95% CIs on every headline metric (resample variants).
- [x] Calibration (ECE + reliability curves) for all classifiers.
- [x] **Negative controls**: shuffle labels → AUROC must collapse to ~0.5.
- [x] **Source-holdout** and **trait-holdout** reported alongside variant-holdout.
- [ ] Multiple-testing correction for pathway enrichment (BH-FDR).

**Agent prompt:**
```
Create scripts/evaluate_suite.py producing results/tables/headline_metrics_ci.csv with
bootstrap 95% CIs (1000 resamples over variant_id) for each feature_group x target x split,
reliability-curve PNGs, a label-permutation negative control (assert AUROC_perm ~ 0.5),
and BH-FDR-adjusted pathway enrichment vs a matched-MAF background. Emit results/figures/
fig_calibration.png and fig_negative_control.png.
```

**Harness** (`tests/gate_eval.py`):
```python
import pandas as pd
ci = pd.read_csv("results/tables/headline_metrics_ci.csv")
neg = ci[ci.feature_group=="PERMUTED"]
assert (neg["AUROC_hi"] >= 0.5).all() and (neg["AUROC_lo"] <= 0.55).all(), "neg control not ~0.5"
real = ci[(ci.feature_group=="AF_LD")&(ci.split=="split_variant")]
assert real["AUROC_lo"].max() > 0.5, "AF_LD CI must exclude chance"
print("GATE PASS: CIs + negative control valid")
```

---

## Phase 6 — Reproducibility & packaging

- [ ] Snakemake DAG wiring all stages (download → features → labels → model → eval → figures).
- [x] `config/` YAMLs (traits, PGS IDs, sources, thresholds) — no magic numbers in code.
- [x] Data registry (`data/registry.yaml`) with URLs, md5s, builds, licenses.
- [x] `pytest` suite = all gate_*.py above + schema validators.
- [ ] `environment.yml` pinned; `README` quickstart updated; DOI via Zenodo.

**Agent prompt:**
```
Create workflow/Snakefile with one rule per stage and correct input/output wiring so
`snakemake -j4 all` reproduces every artifact from scratch. Move hardcoded PGS IDs, trait
maps, source lists, and thresholds into config/atlas.yaml. Add data/registry.yaml recording
each external file's URL, sha256, genome build, and license. Convert all tests/gate_*.py into
a pytest suite and add schema checks (column names/dtypes) for every parquet under data/.
```

**Harness:** `snakemake -n all` (dry-run) resolves with no missing inputs; `pytest -q` green.

---

## Phase 7 — Paper / deliverables

- [ ] Main figures: atlas heatmap, distance-sensitivity, mechanism ablation,
      generalization gap, calibration, pathway-risk, ancestry/pathway coverage gaps (Q4).
- [ ] Results tables with CIs; methods (data provenance, harmonization, splits);
      limitations (no individual phenotypes → behavior not accuracy).
- [ ] Pre-register the 4 research questions Q1–Q4 and map each to a figure/table.

**Agent prompt:**
```
Create paper/figures.py that assembles publication-ready multi-panel figures from
results/tables, and paper/results_tables.py that emits LaTeX/CSV tables with bootstrap CIs.
Draft paper/RESULTS.md mapping Q1->atlas+distance, Q2->pathway-risk, Q3->source-holdout
rank instability, Q4->coverage gaps, each citing the exact artifact path.
```

---

---

## Phase 8 — Genome-wide Scaling (Non-Negotiable)

> **This is the first thing reviewers will attack.** A chr22-only analysis is a
> prototype. Every result below must replicate genome-wide before the paper is ready.

### 8.1 Scale all pipeline stages to chromosomes 1–22

- [ ] chr1–22 genotype preprocessing
- [ ] chr1–22 PGS harmonization
- [ ] chr1–22 score-shift atlas
- [ ] chr1–22 AF features
- [ ] chr1–22 LD features
- [ ] chr1–22 SEL / constraint features
- [ ] chr1–22 GWAS concordance labels (Pan-UKBB + BBJ + FinnGen)
- [ ] chr1–22 pathway mapping

**Agent prompt:**
```
Extend scripts/run_genomewide.py to loop --chrom 1..22 across every feature and label
script, writing per-chrom parquets then concatenating. Keep per-chromosome memory bounded.
Verify >=20 chromosomes appear in master_variant_table_genomewide.parquet and len > 1M rows.
```

**Minimum journal-ready gate:**
- AF+LD+SEL AUROC must replicate genome-wide (≥0.80 on associated variant-holdout).
- Trait-holdout AUROC must remain ≈chance (≤0.55) genome-wide, or the deviation must be explained.
- Pathway enrichments must survive genome-wide LD-aware correction.

### 8.2 New experiments enabled by genome-wide data

Run these only after 8.1 is complete:

- [ ] **Chromosome-holdout:** train on 21 chromosomes, test on held-out chromosome. Report per-chromosome AUROC.
- [ ] **LD-block holdout:** split by ~independent LD blocks (use LDetect block files), not random variants.
- [ ] **Per-chromosome replication:** report AF/LD/SEL AUROC separately for each of 22 chromosomes; flag outliers.
- [ ] **Leave-one-locus-out (LOLO):** confirm results are not driven by a single high-signal locus; re-run ablation excluding the top-contributing locus by SHAP value.
- [ ] **Genome-wide pathway enrichment:** redo Reactome / GO / KEGG analysis across all chromosomes with LD-block permutation and BH-FDR.
- [ ] **Exclude complicated regions:** repeat all main analyses with and without MHC (chr6:28–34 Mb), known large inversions (chr8p23, chr17q21), and long-range LD regions. Report delta-AUROC.

**Resources:**
- LDetect LD blocks: `https://bitbucket.org/nygcresearch/ldetect-data` (GRCh38 per-pop block files)
- MHC coordinates: chr6:28,510,020–33,480,577 (GRCh38)

---

## Phase 9 — Strengthen the Dependent Variable

> Reviewers will ask: "Is your I² measuring true transferability, or just noise,
> sample-size imbalance, and imperfect allele harmonization?"

### 9.1 Add covariates that control for measurement artifacts

- [ ] Add per-variant/source metadata to label table: GWAS N, SE, MAF, MAC, imputation INFO score (if available), case-control ratio, effect-allele harmonization flag, distance to nearest gene, VEP consequence, PGS weight magnitude.
- [ ] Refit the concordance label model conditioning on these covariates; check that I² signal survives.

**Agent prompt:**
```
Extend build_concordance_labels_multisource.py to attach per-source metadata columns:
gwas_n, se, maf, mac, info_score (NaN if unavailable), cc_ratio, harmonization_ok,
dist_to_gene, vep_consequence, pgs_weight. Save to gwas_concordance_labels_multisource.parquet.
Add a diagnostic: for each column, report its Spearman correlation with I2 across associated
variants and write to results/tables/label_covariate_audit.csv.
```

### 9.2 Add field-standard cross-ancestry correlation metrics

- [ ] Run **Popcorn** (transethnic genetic correlation) for each trait × source pair; report `rg_EUR_EAS`, `rg_EUR_AFR`, `rg_EUR_SAS`.
- [ ] Run **S-LDSC** (stratified LD-score regression) to get partitioned heritability and annotation enrichment for each trait.
- [ ] **Validate Atlas against Popcorn:** test whether higher mean Atlas portability-risk predicts lower Popcorn rg. Report Spearman r.
- [ ] **Validate Atlas against S-LDSC:** test whether high-risk variants are enriched in S-LDSC annotations with high heritability contribution.

**Resources:**
- Popcorn: `https://github.com/brielin/Popcorn` (pip installable; needs per-pop sumstats + LD scores)
- S-LDSC: `https://github.com/bulik/ldsc` (pip installable; needs LD scores from `https://data.broadinstitute.org/alkesgroup/LDSCORE/`)
- 1000G LD scores by ancestry: `https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_EAS_baselineLD_v2.2_ldscores.tgz` (and EUR, AFR equivalents)

**New metrics to report:**

| Level | Metric |
|---|---|
| Variant | AUROC / AUPRC for high-I² label; Pearson/Spearman r for continuous I²; ECE / Brier; source-adjusted AUROC |
| Trait | Correlation of Atlas risk with Popcorn rg; correlation with cross-ancestry beta correlation; correlation with observed score shift |
| Pathway | Mean I²; high-risk enrichment OR; pathway FDR q-value; correlation with S-LDSC annotation enrichment |

---

## Phase 10 — Fine-mapping to Separate Tag-SNP Heterogeneity from Causal Heterogeneity

> High I² at a tag SNP may mean (a) the tag fails to capture the causal variant in
> another ancestry (LD/tagging artifact) OR (b) the causal effect itself differs.
> These are biologically distinct. You must distinguish them.

### 10.1 Run SuSiE / FINEMAP per trait × source × ancestry

- [ ] For each trait and GWAS source, define loci around genome-wide significant or PGS-index variants (±500 kb).
- [ ] Build ancestry-specific LD matrices from 1000G matched reference panels.
- [ ] Run **SuSiE** (sum of single effects; fits from z-scores + LD reference) per source/ancestry.
- [ ] Run **FINEMAP** as a cross-check (shotgun stochastic search over causal configurations).
- [ ] Extract credible sets and PIPs; compare cross-ancestry credible-set overlap.

**Resources:**
- SuSiE R package: `https://stephenslab.github.io/susieR/` (also susie_rss for summary stats)
- FINEMAP: `http://www.christianbenner.com/` (v1.4, Linux binary)
- SuSiEx (multi-ancestry fine-mapping): `https://github.com/getian107/SuSiEx`

**Agent prompt:**
```
Create scripts/run_finemapping.py. For each trait and GWAS source: (1) identify loci from
associated variants (clump at r2<0.1, p<5e-8 or PGS-index SNPs), (2) extract z-scores and
1000G-ancestry LD matrices per locus, (3) run susie_rss via rpy2 or subprocess R, (4) extract
credible sets and PIPs, (5) compute cross-ancestry Jaccard overlap between credible sets.
Save results/tables/finemapping_results.parquet [variant_id, trait, source, pip, credible_set_id,
cs_jaccard_EUR_EAS, cs_jaccard_EUR_AFR, causal_signal_shared].
```

### 10.2 New fine-mapping-derived labels

- [ ] `tag_I2`: I² at the PGS/tag SNP (existing).
- [ ] `credible_set_overlap`: Jaccard overlap of credible sets across ancestries.
- [ ] `lead_PIP_discordance`: difference in top-PIP causal candidate across ancestry-specific fine-maps.
- [ ] `causal_signal_shared`: binary — same credible-set signal likely shared vs. not.
- [ ] `tag_vs_causal_failure`: classify each locus as (i) stable tag + stable signal, (ii) unstable tag + stable signal, (iii) stable tag + unstable effect, (iv) unstable tag + unstable signal.

> **Important phrasing note:** Do not claim G×E or epistasis from summary statistics.
> Fine-mapping can help distinguish LD/tagging artifacts from plausible effect heterogeneity,
> but cannot prove environment-specific effects without individual-level phenotype/environment data.
> Phrase results as "credible-set-level heterogeneity consistent with effect instability" —
> not "G×E-driven causal heterogeneity."

---

## Phase 11 — Fix the Trait-Holdout Failure with Functional Annotations

> Trait-holdout AUROC ≈ 0.50 is a strong finding, but it creates a utility problem:
> if the model cannot generalize across traits, what makes it broadly useful?
> The right fix is to give the model biological context — not just a bigger model.

### 11.1 Add functional genomic annotation sources

- [ ] **Ensembl Regulatory Build**: open chromatin, ChIP/ATAC-derived regulatory elements (`https://ftp.ensembl.org/pub/release-110/regulation/homo_sapiens/`)
- [ ] **Roadmap Epigenomics**: 111 reference epigenomes, histone marks, chromatin states (`https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChRCellTypeSpecific/`)
- [ ] **ENCODE cCREs**: candidate cis-regulatory elements (`https://api.wenglab.org/screen_v13/fdownloads/GRCh38-cCREs.bed`)
- [ ] **GTEx v8 eQTL / sQTL annotations**: whether variant is a significant eQTL, in which tissues (`https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression`)
- [ ] **VEP consequence annotations** (already partially available via gnomAD; expand to full VEP output)
- [ ] **CADD v1.7 scores**: `https://cadd.gs.washington.edu/download`
- [ ] **phyloP / phastCons**: see Phase 2.2 above

### 11.2 Trait-holdout experiment with functional annotations

- [ ] Re-run trait-holdout under 5 models:
  - Model A: AF + LD + SEL (current)
  - Model B: AF + LD + SEL + VEP consequence
  - Model C: AF + LD + SEL + Roadmap tissue chromatin states (top 5 trait-relevant tissues)
  - Model D: AF + LD + SEL + Ensembl regulatory features
  - Model E: AF + LD + SEL + GTEx eQTL tissue specificity + gene/pathway embeddings
- [ ] For each model, report trait-holdout AUROC with 95% CI.
- [ ] Key question: does knowledge of vascular / liver / adipose / immune / pancreatic regulatory context recover cross-trait generalization for CAD / LDL / BMI / T2D?

> **Either result is publishable:** if annotations improve trait-holdout, report the recovery.
> If they do not, the conclusion becomes: "portability risk is highly trait-specific and cannot
> be inferred from generic population-genetic or functional features alone."

---

## Phase 12 — Journal-Grade Pathway Analysis

> Current pathway atlas is hypothesis-generating. To be publishable it needs
> LD-aware, gene-density-aware, and annotation-robust pathway tests.

### 12.1 Required pathway controls

- [ ] **BH-FDR correction** across all tested pathways (Reactome + GO + KEGG + MSigDB).
- [ ] **LD-block permutation**: shuffle labels across LD blocks (not individual SNPs) to generate the null.
- [ ] **Matched null pathways**: match each tested pathway by number of SNPs, number of genes, gene length, LD score, MAF, and mean PGS weight; compare observed enrichment to matched null.
- [ ] **Locus-collapsed analysis**: collapse variants into independent loci before pathway testing.
- [ ] **Leave-one-locus-out pathway sensitivity**: for each top pathway, remove the highest-signal locus and recompute; flag pathways that collapse.
- [ ] **Variant→gene mapping robustness**: repeat with nearest gene, ±50 kb, ±100 kb, VEP-assigned gene, and GTEx eQTL-assigned gene.
- [ ] **Annotation database replication**: confirm top pathways replicate in at least two of Reactome, GO BP, KEGG, MSigDB Hallmarks.
- [ ] **PGS-weight robustness**: compare unweighted vs. abs-PGS-weight-weighted pathway risk scores.

### 12.2 Pathway output table requirements

Each pathway must report: mean I², fraction high-risk variants, mean predicted risk, dominant mechanism (AF / LD / SEL), number of independent LD blocks, number of genes, top contributing loci, FDR q-value, leave-one-locus-out stability flag, annotation-database replication count.

> **Strong claim threshold:** only call a pathway "robust" if it survives genome-wide
> analysis, LD-block permutation, BH-FDR, leave-one-locus-out, and replicates in ≥2
> annotation sources. All others are labeled "hypothesis-generating."

---

## Phase 13 — Stronger Baselines

> Reviewers will expect both simple and state-of-the-art comparisons.

### 13.1 Simple univariate baselines

- [ ] Add to ablation: FST only, genetic distance only, MAF only, LD score only, PGS weight magnitude only, GWAS p-value only, GWAS SE only, GWAS N only, chromosome + position only, nearest-gene only, VEP consequence only.

> Purpose: prove the model is not just learning "rare variants have noisy betas" or
> "large-effect variants replicate better."

### 13.2 ML baselines

- [ ] Add: logistic regression, elastic net, random forest, XGBoost / LightGBM, flat MLP; compare all against HistGradientBoosting and HPRN-lite.

### 13.3 PRS method baselines (two comparison modes)

**Mode 1 — predict where PRS-CSx / BridgePRS struggle:**

- [ ] Obtain or run PRS-CSx and BridgePRS weights per ancestry.
- [ ] Define residual weight-instability target: `abs(beta_target_observed - beta_method_predicted)`.
- [ ] Test: does Atlas risk predict residual method error? Report Spearman r.
- [ ] If yes, claim: "The Atlas identifies variants where state-of-the-art multi-ancestry PRS methods still face transferability uncertainty."

**Mode 2 — use Atlas to filter/reweight variants:**

- [ ] Create filtered PGS variants: remove top 5%, 10%, 20% highest-risk variants; compare to random removal and FST-based removal.
- [ ] Create reweighted PGS: `new_weight = old_weight × (1 - predicted_risk)`.
- [ ] Evaluate using **phenotype-free metrics only** (cross-ancestry beta concordance, score distribution shift, top-percentile ancestry composition, credible-set concordance, source-holdout weight error).
- [ ] Do NOT claim clinical improvement without individual-level phenotype data.

**Resources:**
- PRS-CSx: `https://github.com/getian107/PRScsx`
- BridgePRS: `https://github.com/clivehoggart/BridgePRS`

---

## Phase 14 — Source-Holdout as a Core Result

> If your model only works when the same source contributes to both train and test,
> reviewers will say it learned source-specific artifacts. Source-holdout is stronger
> evidence of mechanism learning.

- [ ] **Train: Pan-UKBB + BBJ → Test: FinnGen**
- [ ] **Train: Pan-UKBB + FinnGen → Test: BBJ**
- [ ] **Train: BBJ + FinnGen → Test: Pan-UKBB**
- [ ] **Train: Pan-UKBB EUR/EAS/AFR → Test: external OpenGWAS or FinnGen endpoint** (a trait not used during training)
- [ ] For each: report AUROC, AUPRC, I² Spearman r, calibration ECE, and performance drop vs. random split with 95% CI.

**Resources:**
- OpenGWAS: `https://gwas.mrcieu.ac.uk/` (thousands of publicly available GWAS sumstats)

---

## Phase 15 — Formal Trait-Holdout Experiment

> Make the trait-holdout collapse rigorous, not just a side observation.

- [ ] For each held-out trait (CAD, T2D, BMI, LDL, plus 2+ new traits from Phase 17), train on all remaining traits and test. Report mean AUROC, AUPRC, I² Spearman r, calibration, and permutation control with bootstrap CI.
- [ ] **Trait similarity analysis**: create trait embeddings from pathway enrichment profiles, S-LDSC tissue enrichment, PGS variant annotation distribution, and GWAS genetic correlation. Test: do biologically related traits generalize better?
  - Expected: LDL ↔ CAD should generalize better than BMI ↔ CAD.
  - If no traits generalize: conclude "trait-specific models are necessary" — still publishable.
- [ ] Formally test: is trait-holdout AUROC correlated with cross-trait genetic correlation (rg)?

---

## Phase 16 — Local Ancestry Experiment (Secondary, High-Novelty)

> Not required for the main paper but substantially raises novelty.
> Prior work (Harpak/Wang, Hu et al.) argues global ancestry groupings are insufficient
> and points toward local/trait-specific mechanisms.

- [ ] Use admixed 1000G populations: ASW, ACB, MXL, PEL, PUR, CLM.
- [ ] Run **RFMix2** on chr22 (and genome-wide after Phase 8) to get per-SNP local ancestry calls.
- [ ] For each PGS SNP, extract: global ancestry PCs, local ancestry call, per-SNP PRS contribution, Atlas risk, LD divergence, AF divergence.
- [ ] Test: (a) does local ancestry explain per-locus PRS contribution deviation better than global PCs? (b) do high-Atlas-risk loci show stronger local ancestry dependence? (c) do high-LD-risk variants show stronger local ancestry dependence than high-AF-risk variants?

**Resources:**
- RFMix2: `https://github.com/slowkoni/rfmix`
- 1000G admixed populations: already in the Phase 8 1000G resource above.

---

## Phase 17 — Expand to More Traits (Contrastive Biology)

> Do not add random traits. You need contrastive biology to test cross-trait structure.

- [ ] **Cardiometabolic:** add HDL, triglycerides, fasting glucose (OpenGWAS / Pan-UKBB).
- [ ] **Immune/inflammatory:** add asthma, rheumatoid arthritis, IBD, WBC count, lymphocyte count, eosinophil count (FinnGen + Pan-UKBB both have these).
- [ ] **Anthropometric:** add height, waist-hip ratio.
- [ ] **Neuro (if sufficient cross-ancestry data):** Alzheimer's (IGAP + GCAD), Parkinson's (GP2).

**Key questions enabled by contrastive traits:**
- Do immune traits show stronger selection/turnover signal than metabolic traits?
- Do lipid traits show more conserved causal architecture than immune traits?
- Does trait-holdout improve within biological families (LDL→CAD) vs. across families (T2D→RA)?

---

## Phase 18 — Intervention Experiments

> This is what makes the work **actionable**, not just descriptive.

### 18.1 High-risk variant filtering

- [ ] For each PGS: remove top 5%, 10%, 20% highest-risk variants. Compare to: random removal, highest-FST removal, highest-LD-score removal, lowest-MAF removal.
- [ ] Measure: cross-ancestry beta concordance change, score distribution shift, top-percentile ancestry composition change.

### 18.2 Risk-weighted reweighting

- [ ] `new_weight = old_weight × (1 − predicted_portability_risk)` or exponential decay variant.
- [ ] Measure: score-shift reduction, beta-concordance improvement, fraction of original weight retained, proxy for trait-signal loss.

### 18.3 Per-patient confidence flag

- [ ] For each individual, compute a `score_confidence` = weighted average portability confidence over the variants contributing most to that individual's score (top decile by |weight × genotype|).
- [ ] Report as a reliability flag only. Do not claim clinical calibration without individual-level phenotype data.

**Key claim to establish:** "The Atlas can prospectively identify variants that should be downweighted, excluded, or flagged before cross-ancestry deployment."

---

## Phase 19 — Comprehensive Negative Controls

> These are what prevent every major result from being dismissed.

### Labels
- [ ] Permute I² labels globally
- [ ] Permute I² labels within trait
- [ ] Permute I² labels within chromosome
- [ ] Permute I² labels within LD block
- [ ] Permute sign labels within source pair

### Features
- [ ] Random AF features (same mean/SD, no structure)
- [ ] Random LD features
- [ ] Random selection features
- [ ] Position-only model
- [ ] Chromosome-only model
- [ ] MAF-only model
- [ ] PGS-weight-magnitude-only model
- [ ] SE / GWAS-N-only model

### Pathway
- [ ] Gene-label permutation
- [ ] LD-block-preserving label permutation
- [ ] Matched random gene sets (same size, same chromosome distribution)
- [ ] Pathway-size-matched null
- [ ] SNP-density-matched null
- [ ] Leave-one-gene-out for each top pathway
- [ ] Leave-one-locus-out for each top pathway

---

## Phase 20 — HPRN-lite (Upgraded Spec)

> Supersedes Phase 4. Do not build HPRN until Phase 8 (genome-wide) and Phase 11
> (functional annotations) are complete. The architecture below incorporates all
> lessons from the trait-holdout failure.

### Required architecture components

- [ ] **AF branch**: encodes allele-frequency divergence features per variant
- [ ] **LD branch**: encodes LD divergence and local block structure
- [ ] **SEL / conservation branch**: PBS, LOEUF, phyloP, iHS
- [ ] **Functional annotation branch**: Roadmap chromatin states, GTEx eQTL, VEP consequence, CADD
- [ ] **Gene / pathway aggregation**: Reactome-hierarchical pooling using ReactomePathwaysRelation.txt
- [ ] **Trait embedding**: trait category, tissue relevance, S-LDSC tissue enrichment, GWAS summary architecture
- [ ] **Prediction heads**: `y_high_I2` (BCE), `I2` continuous (MSE), `sign_discordance` (BCE), `causal_signal_shared` (BCE from Phase 10)

### Required ablations

- [ ] Without pathway context
- [ ] Without trait embedding
- [ ] Without functional annotations
- [ ] Without SEL / conservation
- [ ] Without LD branch
- [ ] Without AF branch
- [ ] Flat MLP (same parameter count, no hierarchy)
- [ ] HistGradientBoosting (current best baseline)

### Journal-ready HPRN claim threshold

HPRN is only publishable if it **improves over HistGradientBoosting** on at least three of:
trait-holdout, source-holdout, chromosome-holdout, pathway-level FDR replication, and calibration ECE.
Improvement on random variant-holdout alone is not sufficient.

---

## Reusable harness patterns

- **Smoke test before scale:** every new script takes `--limit N` to run on a small slice
  first; CI runs the limited mode.
- **Gate scripts** (`tests/gate_*.py`): each returns non-zero on failure and prints
  `GATE PASS: ...` on success; chain them in `tests/run_all_gates.py`.
- **Schema validator:** `tests/test_schemas.py` asserts required columns/dtypes per parquet.
- **Leakage check:** assert no `variant_id` appears in >1 fold for `split_variant`.
- **Determinism:** all randomness seeded (`SEED=719`); gates assert identical metrics on rerun.

```bash
# one-shot acceptance run
python tests/run_all_gates.py && python -m pytest tests/ -q
```
