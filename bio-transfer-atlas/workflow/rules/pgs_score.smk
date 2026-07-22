"""
Snakemake rules: apply harmonized PGS weights to 1000G individuals via PLINK2.
"""


rule calculate_scores:
    input:
        pgen       = "data/interim/1000g/chr{chrom}.qc.pgen",
        score_file = "data/interim/pgs/{pgs_id}/chr{chrom}.harmonized.tsv",
    output:
        sscore = "data/processed/scores/{pgs_id}/chr{chrom}.sscore",
    shell:
        """
        mkdir -p data/processed/scores/{wildcards.pgs_id}
        plink2 \
          --pfile data/interim/1000g/chr{wildcards.chrom}.qc \
          --score {input.score_file} 1 2 3 header cols=+scoresums \
          --out data/processed/scores/{wildcards.pgs_id}/chr{wildcards.chrom}
        """


rule merge_scores:
    input:
        expand(
            "data/processed/scores/{{pgs_id}}/chr{chrom}.sscore",
            chrom=config["genotype"]["chromosomes"],
        ),
    output:
        "data/processed/scores/{pgs_id}/all_chrom.sscore",
    script:
        "../../scripts/calculate_scores.py"
