---
license: mit
task_categories:
  - other
language:
  - en
tags:
  - genomics
  - gwas
  - polygenic-scores
  - population-genetics
size_categories:
  - 100G<n<1T
---

# FAIRGEN cold backup layout (Hugging Face)

Repository: [mj0jpgg/fairgen](https://huggingface.co/datasets/mj0jpgg/fairgen)

## What is here

Bulky **regenerable** pipeline artifacts mirrored from `bio-transfer-atlas/data/`:

| Hub path | Local path | ~Size | Regenerate via |
|----------|------------|-------|----------------|
| `backup/raw/` | `data/raw/` | ~47 GB | `scripts/download_*` |
| `backup/interim/` | `data/interim/` | ~138 GB | `scripts/preprocess_grch38_genomewide.py` |
| `backup/features/` | `data/features/` | ~288 GB | `scripts/compute_*_features.py` |

## Kept on disk (not in this backup)

These stay local to resume modeling and scoring without re-download:

- `data/modeling/` — master variant tables, portability model
- `data/processed/` — harmonized PGS weights, score matrices
- `data/labels/` — GWAS concordance labels (~3 GB)
- `data/annotations/` — gene/pathway maps

## Restore

```bash
pip install huggingface_hub
export HF_TOKEN=...  # read access

# Restore one directory
huggingface-cli download mj0jpgg/fairgen --repo-type dataset \
  --include "backup/raw/**" --local-dir ./restore

# Then move restore/backup/raw -> data/raw
```

Or use `upload_hf_backup.py` in reverse with `huggingface_hub.hf_hub_download`.

## Provenance

See `data/registry.yaml` and `data/raw/MANIFEST.tsv` in the GitHub repo for source URLs and licenses.
