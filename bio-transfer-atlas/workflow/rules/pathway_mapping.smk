"""
Snakemake rules: variant -> gene -> Reactome pathway mapping.
"""


rule build_variant_gene_map:
    input:
        score_file = "data/interim/pgs/{pgs_id}/chr{chrom}.harmonized.tsv",
    output:
        v2g = "data/interim/pathway/{pgs_id}/chr{chrom}.v2g.tsv",
    params:
        strategy = config.get("pathway_mapping", {}).get("primary", "reactome_positional_10kb"),
        window_kb = 10,
    shell:
        """
        python -m bta.pathways.v2g \
          --score   {input.score_file} \
          --strategy {params.strategy} \
          --window  {params.window_kb} \
          --out     {output.v2g}
        """


rule map_genes_to_pathways:
    input:
        v2g          = "data/interim/pathway/{pgs_id}/chr{chrom}.v2g.tsv",
        reactome_map = "data/raw/reactome/Ensembl2Reactome.txt",
    output:
        g2p = "data/interim/pathway/{pgs_id}/chr{chrom}.g2p.tsv",
    shell:
        """
        python -m bta.pathways.reactome \
          --v2g     {input.v2g} \
          --reactome {input.reactome_map} \
          --out     {output.g2p}
        """
