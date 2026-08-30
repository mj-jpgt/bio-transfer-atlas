# Production runbook

Run from the `hostbias/` directory on the 30-core Ubuntu VM. The commands are
foreground commands: Snakemake owns locking and restart state; there is no
daemon, PID file, or second scheduler.

## 1. Prepare the shared-storage overlay

The NFS/virtiofs mount must already exist. The command deliberately refuses to
create a missing mount path, use `/`, or place references outside the mount.

```bash
hostbias production-prepare \
  --base-config runtime/config.primary.yaml \
  --nfs-root /data/hostbias \
  --run-id gate-a-20260729 \
  --output-config runtime/operations/gate-a-20260729.yaml \
  --evidence results/aggregate/operations/gate-a-20260729.ready.json
```

This creates `/data/hostbias/scratch/gate-a-20260729/work` with mode `0700`.
The generated overlay remains Git-ignored under `runtime/`; the READY evidence
is aggregate-only and should be committed and pushed. The overlay resolves the
configured relative references to:

```text
/data/hostbias/references/grch38_no_alt.fa
/data/hostbias/references/indexes/grch38_no_alt
/data/hostbias/references/indexes/hprc-balanced-chm13-v1.mmi
/data/hostbias/references/indexes/gtdb-r220-genomes.mmi
```

Do not use `--allow-non-shared-filesystem` in production. That switch exists
only for isolated tests and emergency validation on a known local filesystem.

## 2. Verify the DAG

The launcher refuses a dirty tracked worktree. Commit and push the READY
evidence first, then run:

```bash
hostbias production-launch \
  --config runtime/operations/gate-a-20260729.yaml \
  --stage fetch \
  --dry-run \
  --evidence results/aggregate/operations/gate-a-20260729.fetch.dry-run.json
```

On the frozen 40-sample manifest, clean DAGs contain:

| Target | Cumulative jobs |
|---|---:|
| `fetch_stage` | 42 |
| `normalize_stage` | 82 |
| `filter_stage` | 123 |
| `assemble_stage` | 203 |
| `downstream_bridge` | 606 |

The endpoint target is launched only after the binning lane has produced all 80
`contig_bins.tsv` and `bin_qc.tsv` contract pairs.

## 3. Run stages sequentially

Use the same scheduler settings for every stage:

```bash
hostbias production-launch \
  --config runtime/operations/gate-a-20260729.yaml \
  --stage fetch \
  --cores 30 --jobs 8 --mem-mb 204800 --disk-mb 2500000 \
  --latency-wait-seconds 120 \
  --evidence results/aggregate/operations/gate-a-20260729.fetch.json
```

After a successful status check, repeat with `--stage normalize`, then
`filter`, `assemble`, and `downstream`, using a stage-specific evidence name.
After binning contracts exist, run `--stage endpoint`.

The caps fit the measured VM:

- 30 total cores and at most 8 simultaneous jobs.
- 204,800 MB aggregate declared memory.
- 2,500,000 MB aggregate declared scratch.
- Assembly and competitive mapping use at most 24 cores and 64,000 MB each.
- Host filtering uses at most 16 cores and 32,000 MB.
- MEGAHIT receives its declared 64,000 MB byte limit, not 90% of host RAM.

## 4. Status and exact restart behavior

Status contains counts only—never paths, sample accessions, URLs, checksums, or
sequence data:

```bash
hostbias production-status \
  --config runtime/operations/gate-a-20260729.yaml \
  --output results/aggregate/operations/gate-a-20260729.status.json
```

If SSH disconnects, a tool fails, or the VM restarts, rerun the exact
`production-launch` command for that stage. Every invocation uses:

```text
--rerun-incomplete --keep-going
--rerun-triggers mtime input params code software-env
```

It never uses `--force`, deletes completed outputs, or disables Snakemake's
lock. Complete files are retained, incomplete rule outputs are regenerated, and
independent samples continue after a sample-level failure. Launch evidence is
written as `RUNNING` before Snakemake starts and finalized as `COMPLETE`,
`FAILED`, or `DRY_RUN_COMPLETE`, with exit code and aggregate stage counts.
