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
# Generate runtime/config.sentinel.yaml using docs/DATA_PROVENANCE_RUNBOOK.md.
hostbias validate --config runtime/config.sentinel.yaml
hostbias provenance --config runtime/config.sentinel.yaml \
  --output results/aggregate/provenance.json
snakemake --profile profiles/vm \
  --configfile runtime/config.sentinel.yaml \
  --dry-run all
```

The workflow expects an ENA-resolved runtime manifest and explicit reference
paths. It never downloads or commits controlled-access data. See
`config/config.example.yaml` for the configuration interface and
`docs/DATA_PROVENANCE_RUNBOOK.md` for exact runtime-generation commands.

The six-run dataset-eligibility sentinel is executable independently before the
full workflow. See `docs/SENTINEL_RUNBOOK.md` for the exact restartable VM
command and its aggregate-only output contract.

## Workflow stages

1. Fetch checksum-pinned public paired FASTQs.
2. Trim, synchronize, and deterministically normalize read pairs.
3. Apply GRCh38 host filtering with source-style and strict pair semantics.
4. Assemble retained reads with MEGAHIT.
5. Produce identifier-free assembly QC.
6. Align contigs independently to checksum-verified human and GTDB minimap2
   indexes, then convert PAF records into the competitive-label TSV contract.
7. Combine externally produced bin assignments/QC with competitive labels and
   publish identifier-free sample endpoint aggregates.

Paired-read stages validate checksums or synchronization before publication.
Every compute rule declares threads and memory and is restartable. `all` ends at
assemblies. `downstream_bridge` runs assembly QC and competitive mapping.
`gate_a_endpoint_aggregates` additionally expects standardized binning contracts
under `work/binning/{sample}/{mode}/` and publishes sample aggregates.

```bash
snakemake --profile profiles/vm downstream_bridge
snakemake --profile profiles/vm gate_a_endpoint_aggregates
```
