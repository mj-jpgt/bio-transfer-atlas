"""
Snakemake rules: export processed outputs to data/processed/ for HuggingFace upload.
"""


rule export_processed:
    input:
        "results/tables/population_shift_metrics.parquet",
        "data/processed/ancestry/sample_metadata.parquet",
        "data/processed/ancestry/ancestry_pcs.parquet",
    output:
        touch("results/.export_done"),
    shell:
        """
        echo "Processed outputs ready in data/processed/ and results/tables/"
        """
