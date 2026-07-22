"""
Snakemake rules: compute all instability metrics.
"""


rule compute_population_shift:
    input:
        score_matrix = "data/processed/scores/{pgs_id}/all_chrom.sscore",
        sample_meta  = "data/processed/ancestry/sample_metadata.parquet",
    output:
        shifts = "data/processed/metrics/{pgs_id}/population_shift.parquet",
    script:
        "../../scripts/compute_metrics.py"


rule aggregate_metrics:
    input:
        shifts = expand(
            "data/processed/metrics/{pgs_id}/population_shift.parquet",
            pgs_id=[]  # filled dynamically from scores.yaml at runtime
        ),
        sample_meta = "data/processed/ancestry/sample_metadata.parquet",
    output:
        "results/tables/population_shift_metrics.parquet",
    script:
        "../../scripts/compute_metrics.py"
