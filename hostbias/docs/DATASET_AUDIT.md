# Stage 0 dataset audit

This document freezes the data choices for the initial host-filtration
propagation experiment. It must be read together with
`config/stage0_selection.yaml`; a sample may not be excluded or replaced on the
basis of host rate or any downstream result.

## Primary comparison

| Arm | BioProject | Cohort label used in Stage 0 | Role |
| --- | --- | --- | --- |
| Tanzania | PRJNA686265 | Tanzania cohort | Primary |
| Netherlands | PRJNA319574 | Netherlands cohort | Primary |

Both projects contain open raw paired-end Illumina shotgun metagenomes. Stage 0
estimates a **cohort-origin association**, not a causal ancestry effect. The two
arms differ in geography, collection, laboratory processing, and study design.
Constructed HPRC ground truth is required before making a causal ancestry claim.

The exact 20 primary and 10 reserve runs per arm are frozen in
`config/stage0_samples.tsv`. Runs were ranked within BioProject by ENA
`base_count` descending, with accession ascending as the deterministic
tie-breaker, after applying the pre-result technical filters. A reserve can
replace a primary only for an allowed technical failure.

## Sources rejected for the primary comparison

| Source | Stage 0 disposition | Reason |
| --- | --- | --- |
| saMBA | Excluded | The relevant release is 16S/DADA2 rather than raw shotgun metagenomics. |
| AWI-Gen 2 | Excluded from initial leakage test | Human-containing reads require controlled access; the open split cannot measure initial GRCh38 escape. |
| PRJNA678454 | Replication only | Reads were already filtered against hg19, so this source can test downstream residue but not initial leakage. |

## Frozen fallback order

Fallbacks are activated only if a whole primary source fails the sentinel
eligibility rule:

1. Replace an ineligible Netherlands source with Danish MetaHIT `PRJEB4336`.
2. If the Tanzania source is ineligible, compare 20 Peru and 20 USA samples from
   `PRJNA268964`.
3. Use `PRJNA278393` as an 11-versus-11 Hadza/Italian protocol-matched
   replication.
4. Use `PRJNA678454` only to test propagation after prior hg19 filtering.

Activating a fallback requires a committed eligibility report identifying the
failed criterion. The fallback's run manifest must be frozen before downloading
the full cohort.

## Privacy and repository boundary

Raw reads, reference FASTAs, alignments, read-level calls, and human-derived
contigs remain on restricted VM storage and never enter Git. Only accessions,
public metadata, checksums, aggregate results, figures, and sanitized logs are
versioned. Raw FASTQ files are deleted only after filtration outputs and their
checksums have been verified.

