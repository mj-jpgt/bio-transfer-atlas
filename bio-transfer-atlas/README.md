# Biological Transferability Atlas

Open-data atlas for **cross-ancestry PRS/GWAS portability risk**: predict variant-level failure from AF/LD/selection features under LD-block CV, then test Catalog score edits on 1000 Genomes with honest nulls (matched Monte Carlo, Duffy, PAGE).

**Start here:** the repository root [`../README.md`](../README.md) (full overview, claims, layout).  
**Manuscript:** [`paper/`](paper/) · **Freeze tables:** [`results/tables/`](results/tables/) · [`BTA_robustness_freeze_results_20260720/`](BTA_robustness_freeze_results_20260720/)  
**Gate:** `python scripts/gate_literature_roadmap.py`

## Local quickstart

```bash
mamba env create -f environment.yml
mamba activate bta
make test
make smoke   # optional chr22 Snakemake smoke
```

## Not in this clone

Raw genotypes, Pan-UKB dumps, and large parquet/pgen artifacts are gitignored. Rebuild via `scripts/download_*.py` / Lambda workflows documented in `LOCAL_COMPUTE.md` and `ROADMAP_TO_COMPLETION.md`.

## What this is / is not

| Is | Is not |
|----|--------|
| Phenotype-free concordance risk models | Clinical risk calculator |
| Ancestry **mean separation** (MAD) for edits | Proof of repaired phenotype \(R^2\) |
| Descriptive GAT / signed SuSiE tiers | Mechanism “solved by” graphs |

Primary AUROC citations must use **LD-block** splits (≈0.627 for AF_LD_SEL).
