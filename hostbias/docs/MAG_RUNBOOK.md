# MAG binning, QC, and taxonomy runbook

## Design

The primary Gate A MAG result is the non-redundant DAS Tool consensus of three
independent binners:

- MetaBAT2 2.17
- MaxBin2 2.2.7 with the bacteria/archaea-aware 40-marker set
- CONCOCT 1.1.0

All three receive the same contigs eligible at the frozen 1,500 bp threshold
and coverage derived by mapping that sample's retained read pairs back to its
own assembly. This avoids introducing cohort-dependent co-abundance panels.
The three original bin maps and bin directories remain available under
`work/binning/` for method-sensitivity analysis; DAS Tool 1.1.7 is the
preregistered primary endpoint.

Consensus bins are assessed without pre-filtering by CheckM2 1.1.0 and GUNC
1.1.1, then classified with GTDB-Tk 2.4.0 and the frozen R220 database.
The three tools use separate pinned environments because current CheckM2 and
GUNC require incompatible DIAMOND versions.
Filtering occurs only in the endpoint definition. A GUNC result of `nan` is
translated to `gunc_pass=false`. Missing tool rows are fatal. A valid consensus
with zero bins produces header-only private contracts and a zero-propagation
sample endpoint.

## Database variables

Set these variables to readable, immutable database installations. They
override the paths in the runtime config:

```bash
export CHECKM2DB=/data/hostbias/databases/checkm2/uniref100.KO.1.dmnd
export GUNC_DB=/data/hostbias/databases/gunc/gunc_db_progenomes2.1.dmnd
export GTDBTK_DATA_PATH=/data/hostbias/databases/gtdbtk/release220
```

`CHECKM2DB` and `GUNC_DB` must point to DIAMOND database files, not their parent
directories. `GTDBTK_DATA_PATH` must point to the unpacked GTDB-Tk R220 data
root. GTDB-Tk 2.4.0 is pinned because it requires R220. The shared preflight
validates that all three database payloads are present; each tool then performs
its own format/install validation inside its dedicated environment.

## Live targets

From the `hostbias` project directory:

```bash
snakemake --profile profiles/vm \
  --configfile runtime/config.primary.yaml \
  mag_database_preflight

snakemake --profile profiles/vm \
  --configfile runtime/config.primary.yaml \
  mag_consensus

snakemake --profile profiles/vm \
  --configfile runtime/config.primary.yaml \
  mag_endpoint_inputs

snakemake --profile profiles/vm \
  --configfile runtime/config.primary.yaml \
  gate_a_endpoint_aggregates
```

The targets are restartable in that order, but the final target can also be run
directly and Snakemake will materialize all prerequisites.

## Private/public boundary

All BAMs, depth tables, raw binner outputs, DAS Tool bins, CheckM2/GUNC/GTDB-Tk
reports, and standardized `contig_bins.tsv`/`bin_qc.tsv` files stay beneath
`work/`. No MAG sequence or contig/bin identifier is written beneath
`results/`. Only the existing identifier-free sample endpoint JSON is
publishable.
