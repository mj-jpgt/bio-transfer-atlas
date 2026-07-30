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
```

On 2026-07-30, the canonical ENA responses contained 320 rows for
`PRJNA686265` and 875 rows for `PRJNA319574`. Technical filtering and
one-run-per-BioSample selection retained 320 and 471 eligible rows,
respectively. Their hashes are frozen in `config/ena_snapshot_ledger.tsv`.
Rebuilding from those snapshots produced a byte-identical 60-run manifest with
SHA-256
`277a2d2633a3fb563510d2a00ee8c33dd0b9378fce30df10ab24d92753dc4c0d`.

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

