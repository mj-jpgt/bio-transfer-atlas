# Local compute budget (16 GB RAM / Windows / OneDrive)

Hard rules for this machine:
1. NEVER download a full Pan-UKB chromosome bgz (chr1 alone ~2 GB compressed). Use scripts/download_panukbb_region.py.
2. NEVER load master_variant_table_genomewide_genomewide.parquet into pandas. Stream-sample to disk (<=500k rows) via scripts/_stream_sample_associated.py.
3. Cap HistGB train at <=250k rows; max_bins=64; early_stopping on.
4. Run ONE heavy job at a time. Delete *.bgz / _tmp* downloads after parquet write.
5. Prefer cached samples under data/modeling/_tmp_*_sample.parquet for re-evals.
6. External PAGE/GBMI: use associations API or chr-filtered streams, not 1 GB FTP dumps.
7. If avail RAM < 3 GB, stop and free disk/OneDrive stubs before continuing.

Primary paper numbers already on disk (do not recompute unless needed):
- LD-block AF_LD_SEL AUROC 0.627
- MHC sensitivity table
- finemap_tiers_genomewide_zlead.parquet (unconfounded)
