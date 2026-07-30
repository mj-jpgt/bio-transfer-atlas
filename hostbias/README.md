# Hostbias Gate A

Hostbias tests whether human reads that evade a conventional GRCh38 host filter
can assemble and propagate into apparently novel microbial genome bins. Gate A
compares identically processed public shotgun metagenomes from Tanzania and the
Netherlands. It measures a cohort association; it does not treat cohort as a
proxy for individual genetic ancestry.

This directory contains a reproducible Snakemake workflow and a small Python
package for validation, provenance, and deterministic result generation.

## Safety boundary

Raw reads, references, alignments, assemblies, bins, and other sequence-bearing
files are local-only. They live under `data/`, `references/`, `resources/`,
`work/`, or non-aggregate portions of `results/`, all of which are ignored by
Git. Only schemas, manifests without sequence data, logs stripped of read-level
content, and aggregate tables/reports under `results/aggregate/` may be
committed.

## Quick start

```bash
cd hostbias
mamba env create -f envs/control.yaml
mamba activate hostbias-control
python -m pip install -e .
hostbias validate --config config/config.yaml
hostbias provenance --config config/config.yaml --output results/aggregate/provenance.json
snakemake --profile profiles/vm --dry-run gate_a
```

The workflow expects an explicit sample manifest and reference paths. It never
downloads or commits controlled-access data. See `config/config.example.yaml`
for the complete interface.

## Workflow stages

1. Fetch checksum-pinned public paired FASTQs.
2. Trim, synchronize, and deterministically normalize read pairs.
3. Apply GRCh38 host filtering with source-style and strict pair semantics.
4. Assemble retained reads with MEGAHIT.
5. Hand assemblies to the separately versioned labelling/QC stages.

Every rule writes atomically, declares threads and memory, and is restartable.
The `gate_a` target is intentionally blocked until the downstream labelling,
binning, statistics, and verdict slices are present.
