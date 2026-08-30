# Competitive human reference panel

Gate A labels residual assembled contigs against `hprc-balanced-chm13-v1`, a
linear competitive panel containing verified T2T-CHM13v2.0 plus two haplotypes
from one HPRC donor in each 1000 Genomes superpopulation. A second donor in
every group is downloaded but held out of the union for control evaluation.

## Frozen population design

The selection never infers population from filenames. It joins, by exact donor
ID:

1. The HPRC Release 2 [assembly index](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/5a939042026331a823a6307fe36a3d7e0188a6e0/data_tables/assemblies_release2_v1.0.index.csv).
2. The HPRC Release 2 [sample metadata](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/5a939042026331a823a6307fe36a3d7e0188a6e0/data_tables/sample/hprc_release2_sample_metadata.csv).
3. The IGSR/1000 Genomes [2013 sample panel](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel).

Eligibility requires matching HPRC and IGSR population codes, an IGSR
superpopulation in `{AFR, AMR, EAS, EUR, SAS}`, and exactly two released HPRC
haplotypes with GenBank accessions and checksum URLs. Donors are sorted by ID
within superpopulation; rank 1 is primary and rank 2 is held out.

| Group | Primary donor (population; haplotype accessions) | Held-out donor (population; haplotype accessions) |
| --- | --- | --- |
| AFR | HG02583 (GWD; GCA_042076985.1, GCA_042076715.1) | HG02922 (ESN; GCA_042031165.1, GCA_042030255.1) |
| AMR | HG01167 (PUR; GCA_044167135.1, GCA_044166785.1) | NA19682 (MXL; GCA_044165125.1, GCA_044164665.1) |
| EAS | HG02040 (KHV; GCA_042032555.1, GCA_042031805.1) | HG02155 (CDX; GCA_043304525.1, GCA_043304655.1) |
| EUR | HG00097 (GBR; GCA_044165215.1, GCA_044164745.1) | HG00099 (GBR; GCA_042032525.1, GCA_042031945.1) |
| SAS | HG03742 (ITU; GCA_042034325.1, GCA_042034195.1) | HG03784 (ITU; GCA_044166095.1, GCA_044165595.1) |

HPRC therefore has usable donor-level public labels when joined to IGSR; the
preregistered alternative-panel fallback is not activated.

## Frozen sizes and integrity

`config/competitive_human_panel.tsv` records every URL, authoritative MD5 URL,
compressed byte count, GenBank accession, and role:

- T2T-CHM13v2.0: 932,696,125 bytes; MD5
  `9280657210e4161147cbe13b022225b9`; SHA-256
  `a2f8b712a1958a41b590ce67743c324d64462e17b6b4fdc62747ff4de4b09c75`.
- Ten primary HPRC haplotypes: 8,925,223,631 compressed bytes.
- Ten held-out HPRC haplotypes: 8,988,374,953 compressed bytes.
- Complete download: 21 files and 18,846,294,709 bytes.
- Union inputs: 11 files and 9,857,919,756 compressed bytes.

The acquisition command re-downloads a file unless its expected size and MD5
match. CHM13 must additionally match its frozen SHA-256. Observed SHA-256
values for every HPRC assembly are recorded in the aggregate checkpoint.
Partial HTTP downloads use Range requests and remain under the ignored
reference root until they can be atomically published.

## Exact VM command

From the `hostbias/` project directory:

```bash
micromamba env update -n hostbias-downstream -f envs/downstream.yaml
mkdir -p /data/hostbias/references
chmod 700 /data/hostbias/references

micromamba run -n hostbias-downstream hostbias reference-build \
  --metadata-sources config/reference_metadata_sources.tsv \
  --donors config/hprc_balanced_donors.tsv \
  --panel config/competitive_human_panel.tsv \
  --reference-root /data/hostbias/references \
  --checkpoint results/aggregate/checkpoints/P15_competitive_human_reference_panel.json \
  --threads 32 \
  --index-batch 64G
```

The command:

1. Downloads and hash-verifies the three pinned metadata files.
2. Reproduces the exact ten-donor selection and checks every HPRC assembly row.
3. Downloads or reuses all 21 verified compressed references.
4. Builds a union from CHM13 plus the five primary diploid donors, prefixing
   every FASTA header with its frozen `reference_id`.
5. Runs:

   ```text
   minimap2 -x asm5 -I 64G -t 32 -d hprc-balanced-chm13-v1.mmi.partial hprc-balanced-chm13-v1.fa
   ```

6. Atomically publishes the index and writes aggregate-only evidence.

The resulting paths are:

```text
/data/hostbias/references/panel/hprc-balanced-chm13-v1.fa
/data/hostbias/references/indexes/hprc-balanced-chm13-v1.mmi
results/aggregate/checkpoints/P15_competitive_human_reference_panel.json
```

The committed checkpoint contains metadata hashes and authorities, exact donor
selection, accessions, compressed sizes, expected/observed checksums, union
sequence/base counts and SHA-256, minimap2 version/parameters, and index
size/SHA-256. It contains no absolute VM paths or reference sequence.

