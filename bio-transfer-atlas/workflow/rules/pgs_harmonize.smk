"""
Snakemake rules: harmonize PGS scoring files against 1000G build.
"""


rule harmonize_pgs:
    input:
        score_file = "data/raw/pgs_catalog/scores/{pgs_id}/{pgs_id}_hmPOS_GRCh37.txt.gz",
        pvar       = "data/interim/1000g/chr{chrom}.qc.pvar",
    output:
        harmonized = "data/interim/pgs/{pgs_id}/chr{chrom}.harmonized.tsv",
        report     = "data/interim/pgs/{pgs_id}/chr{chrom}.match_report.tsv",
    shell:
        """
        python -m bta.pgs.harmonize \
          --score {input.score_file} \
          --pvar  {input.pvar} \
          --out   {output.harmonized} \
          --report {output.report}
        """
