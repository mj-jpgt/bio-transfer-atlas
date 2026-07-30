# Stage 0 sentinel runbook

This command tests the first three frozen primary runs in each arm before any
full-cohort download. It requests SRA spots 1 through 1,000,000 for each run,
validates that fastq-dump produced exactly two synchronized mate files, and
streams Bowtie2 SAM output directly into an aggregate mapped-pair counter.
No SAM or BAM is created.

The legacy `fastq-dump` command is intentional here: SRA Toolkit 3.1.1
supports the preregistered `-N`/`-X` spot interval only through `fastq-dump`;
`fasterq-dump` cannot cap a remote run to the first 1,000,000 spots.

## One-time VM preparation

Run from the `hostbias/` project directory:

```bash
micromamba env create -f envs/preprocess.yaml
mkdir -p /data/hostbias/scratch/sentinel
chmod 700 /data/hostbias/scratch/sentinel
```

The configured GRCh38 Bowtie2 prefix must already have a passing reference
checkpoint. The runner hashes every index shard and includes their combined
SHA-256 in each checkpoint.

## Exact live invocation

```bash
cd /data/hostbias/repo/hostbias
micromamba run -n hostbias-preprocess hostbias sentinel-run \
  --manifest config/stage0_samples.tsv \
  --thresholds config/thresholds.yaml \
  --grch38-index /data/hostbias/references/indexes/grch38_no_alt \
  --scratch-root /data/hostbias/scratch/sentinel \
  --output-dir results/aggregate/sentinel \
  --threads 16
```

The command always selects primary ranks 1, 2, and 3 from each arm in the
frozen manifest. Accessions cannot be supplied on the command line.

## Restart and result semantics

Each successful run creates
`results/aggregate/sentinel/runs/<accession>.json`. On restart, a checkpoint is
reused only when its accession, arm, rank, requested spot count, manifest hash,
tool versions, and GRCh38 index hash still match. Stale or malformed
checkpoints are rerun.

The final `sentinel_eligibility.json` applies the preregistered rule:

- at least 100 GRCh38-mapped pairs per million observed pairs;
- at least two of three runs passing in each arm;
- metadata and before/after FASTQ checksums valid for every run.

Exit code 0 means all six runs completed. Eligibility can still be `false`,
which is a valid scientific dataset-screening result. Exit code 2 means an
operational failure; fix the recorded stage and rerun the same command.

On success and failure, the per-accession scratch directory is deleted.
Successful JSON contains only accession/arm/rank, counts, rates, timings, tool
versions, reference and ephemeral-input checksums, and privacy assertions.
Failure diagnostics contain only a stage, return code, stderr hash/size, and a
sanitized short stderr tail. Sequence-like strings are redacted. Raw reads,
SAM, and BAM are never written beneath the aggregate result directory.

