# Methods/statistics interface

`hostbias analyze` consumes aggregate, non-sensitive TSVs and produces the full
methods-stage audit trail. Raw reads, BAMs, references, and human sequences are
outside this interface and must never be copied into the output directory.

Required input tables have exact headers demonstrated in `tests/fixtures/`:

- `alignments.tsv`: best human/GTDB alignment candidates for each contig.
- `bins.tsv`: at most one selected bin assignment per sample/contig.
- `bin_qc.tsv`: one DAS Tool, CheckM2, GUNC, and GTDB-Tk record per bin.
- `sample_groups.tsv`: one frozen cohort label per expected sample.
- `control_alignments.tsv` and `control_truth.tsv`: held-out HPRC/GTDB controls.
- `sensitivities.tsv`: aggregate rerun means for strict pair and identity analyses.

The threshold file has five optional top-level mappings: `labeling`, `endpoints`,
`controls`, `statistics`, and `gate`. Unknown names are errors. Omitted values use
the preregistered defaults. CheckM2 completeness/contamination, identity,
coverage, rates, and fractions are represented on the `[0, 1]` scale.
Every assembled contig must have a row. Contigs with no human or GTDB hit use
`target_domain=none` and zero for all alignment metrics.

Example:

```bash
hostbias analyze \
  --alignments results/alignments.tsv \
  --contig-bins results/bins.tsv \
  --bin-qc results/bin_qc.tsv \
  --sample-groups config/stage0_groups.tsv \
  --control-alignments results/control_alignments.tsv \
  --control-truth config/control_truth.tsv \
  --sensitivities results/sensitivities.tsv \
  --thresholds config/thresholds.yaml \
  --output-dir results/gate_a
```

The command always writes intermediate calls, sample endpoints, bin human
fractions, controls, statistics, effective thresholds, and sibling JSON/Markdown
verdicts. `OPERATIONAL_FAILURE` means controls, sample counts, or sensitivity
inputs are incomplete; only valid analyses receive a scientific `PASS` or `FAIL`.

## Assembly-to-analysis bridge

The restartable `downstream_bridge` Snakemake target runs minimap2 `asm5`
against independently checksum-verified balanced-human and GTDB indexes. PAFs,
run specs, contig IDs, and the resulting `alignments.tsv` remain under `work/`.
Only aggregate assembly QC and mapping manifests are written below
`results/aggregate/`; their JSON Schemas explicitly prohibit paths and
sequence-derived identifiers.

The `mag_endpoint_inputs` target supplies exact `ContigBinRow` and `BinQcRow`
TSVs at `work/binning/{sample}/{mode}/contig_bins.tsv` and `bin_qc.tsv`.
The `gate_a_endpoint_aggregates` target then publishes privacy-safe sample
endpoints. See `docs/MAG_RUNBOOK.md` for database variables and live commands.

The equivalent commands for one unit are:

```bash
hostbias assembly-qc --assembly final.contigs.fa --sample-id T01 \
  --filter-mode source --output assembly_qc.json
hostbias build-alignment-table --spec run_spec.json \
  --output-tsv alignments.tsv --output-manifest alignment_manifest.json
hostbias aggregate-endpoint --alignments alignments.tsv \
  --contig-bins contig_bins.tsv --bin-qc bin_qc.tsv --sample-id T01 \
  --filter-mode source --thresholds config/thresholds.yaml \
  --output endpoint.json
```
