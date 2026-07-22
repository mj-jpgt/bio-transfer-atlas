"""
Snakemake rules: build ancestry representations from PCA + 1000G labels.
"""


rule build_sample_metadata:
    input:
        panel    = "data/raw/1000g/metadata/integrated_call_samples_v3.20130502.ALL.panel",
        eigenvec = expand(
            "data/processed/1000g/chr{chrom}_pca.eigenvec",
            chrom=config["genotype"]["chromosomes"],
        ),
    output:
        sample_meta = "data/processed/ancestry/sample_metadata.parquet",
        ancestry_pcs = "data/processed/ancestry/ancestry_pcs.parquet",
    script:
        "../../scripts/compute_pcs.py"
