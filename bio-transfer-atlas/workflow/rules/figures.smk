"""
Snakemake rules: generate figures.
"""


rule atlas_heatmap:
    input:
        "results/tables/population_shift_metrics.parquet",
    output:
        "results/figures/atlas_heatmap.png",
    script:
        "../../scripts/make_figures.py"
