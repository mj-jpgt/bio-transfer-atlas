# Data provenance runbook

Run these commands from the `hostbias/` directory with the package installed,
or set `PYTHONPATH=src`. The snapshots belong in the run's provenance directory;
the read URLs they contain are public, but the large read files themselves
belong only on restricted scratch storage.

```bash
python -m hostbias.data_manifest fetch-ena \
  --project PRJNA686265 \
  --output provenance/ena/PRJNA686265.ena.tsv

python -m hostbias.data_manifest fetch-ena \
  --project PRJNA319574 \
  --output provenance/ena/PRJNA319574.ena.tsv

python -m hostbias.data_manifest select \
  --arm tanzania=PRJNA686265=provenance/ena/PRJNA686265.ena.tsv \
  --arm netherlands=PRJNA319574=provenance/ena/PRJNA319574.ena.tsv \
  --output provenance/stage0_samples.rebuilt.tsv

python -m hostbias.data_manifest validate \
  --manifest config/stage0_samples.tsv \
  --report provenance/stage0_manifest_validation.json

python -m hostbias.metadata_audit \
  --manifest config/stage0_samples.tsv \
  --arm tanzania=PRJNA686265 \
  --arm netherlands=PRJNA319574 \
  --snapshot-dir provenance/ena \
  --output results/evidence/live_ena_audit.json
```

On 2026-07-30, the canonical ENA responses contained 320 rows for
`PRJNA686265` and 875 rows for `PRJNA319574`. Technical filtering and
one-run-per-BioSample selection retained 320 and 471 eligible rows,
respectively. Their hashes are frozen in `config/ena_snapshot_ledger.tsv`.
Rebuilding from those snapshots produced a byte-identical 60-run manifest with
SHA-256
`277a2d2633a3fb563510d2a00ee8c33dd0b9378fce30df10ab24d92753dc4c0d`.
The live audit exits nonzero if current ENA metadata no longer reproduces the
freeze. Its aggregate/hash-only JSON is safe to commit; the full snapshots stay
in the run provenance directory.

## Generate executable runtime inputs

The selection snapshots above intentionally reproduce the frozen ranking hash.
Execution needs one additional ENA field, `fastq_bytes`, so fetch canonical
runtime snapshots separately:

```bash
python -m hostbias.data_manifest fetch-runtime-ena \
  --project PRJNA686265 \
  --output provenance/ena/PRJNA686265.runtime.ena.tsv

python -m hostbias.data_manifest fetch-runtime-ena \
  --project PRJNA319574 \
  --output provenance/ena/PRJNA319574.runtime.ena.tsv
```

Join those snapshots to the frozen selection. Start with the six-run sentinel;
generate the 40-run primary scope with the same inputs only after the sentinel
eligibility decision.

```bash
hostbias prepare-runtime \
  --snapshot tanzania=provenance/ena/PRJNA686265.runtime.ena.tsv \
  --snapshot netherlands=provenance/ena/PRJNA319574.runtime.ena.tsv \
  --scope sentinel \
  --evidence results/aggregate/checkpoints/P13_runtime_manifest_sentinel.json

hostbias validate --config runtime/config.sentinel.yaml
snakemake --profile profiles/vm \
  --configfile runtime/config.sentinel.yaml \
  --dry-run all
```

For the primary scope, replace `sentinel` with `primary` and write evidence to
`P14_runtime_manifest_primary.json`. The generated `runtime/` files and source
snapshots under `provenance/` contain public FASTQ locations and checksums, so
both directories are Git-ignored. The evidence JSON contains only ordered run
accessions, aggregate byte counts, and content hashes and is safe to commit.

Runtime generation fails if an accession is missing or duplicated, snapshot
order is non-canonical, frozen metadata has drifted, a run does not have exactly
two FASTQs/checksums/sizes, or a URL contains credentials or query tokens.
Downloads are accepted only when both ENA byte counts and MD5 checksums match.

The 2026-07-30 live runtime snapshots resolved 6 sentinel runs to
26,543,719,200 compressed FASTQ bytes and 40 primary runs to
153,851,026,418 compressed FASTQ bytes. Their exact hashes and ordered
accessions are recorded in the two checkpoint JSON files.

For the sentinel check, create a TSV with these columns:

```text
run_accession arm streamed_spots grch38_mapped_pairs metadata_ok checksum_ok
```

Then run:

```bash
python -m hostbias.sentinel \
  --manifest config/stage0_samples.tsv \
  --metrics provenance/sentinel_metrics.tsv \
  --output results/sentinel_eligibility.json
```

For the human-labelling panel, first construct an HPRC catalog with
`donor_id`, `population_group`, `assembly_url`, and an authoritative `sha256`.
The builder selects donors by group and donor ID, takes equal primary counts,
and holds out the next donor in every group:

```bash
python -m hostbias.reference_manifest build \
  --hprc-catalog provenance/hprc_catalog.tsv \
  --per-group 4 \
  --holdout-per-group 1 \
  --chm13-url "$CHM13_URL" \
  --chm13-sha256 "$CHM13_SHA256" \
  --output provenance/reference_manifest.tsv

python -m hostbias.reference_manifest verify \
  --manifest provenance/reference_manifest.tsv \
  --reference-root /data/hostbias/refs \
  --report provenance/reference_checksum_report.json
```

Do not proceed to cohort labelling until reference verification reports
`"valid": true`. The exact HPRC group mapping and checksums must come from the
frozen HPRC catalog; they must never be inferred from cohort outcomes.
