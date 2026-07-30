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
