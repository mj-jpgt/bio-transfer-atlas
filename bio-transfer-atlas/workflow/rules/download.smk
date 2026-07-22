"""
Snakemake download rules for bio-transfer-atlas.
"""

BASE_1000G = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
VCF_TEMPLATE = "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"


rule download_1000g_metadata:
    output:
        "data/raw/1000g/metadata/integrated_call_samples_v3.20130502.ALL.panel"
    shell:
        """
        mkdir -p data/raw/1000g/metadata
        wget -q -O {output} \
          {BASE_1000G}/integrated_call_samples_v3.20130502.ALL.panel
        """


rule download_1000g_vcf:
    output:
        vcf = "data/raw/1000g/vcf/" + VCF_TEMPLATE,
        tbi = "data/raw/1000g/vcf/" + VCF_TEMPLATE + ".tbi",
    params:
        vcf_url = BASE_1000G + "/" + VCF_TEMPLATE,
    shell:
        """
        mkdir -p data/raw/1000g/vcf
        wget -q -O {output.vcf} {params.vcf_url}
        wget -q -O {output.tbi} {params.vcf_url}.tbi
        """


rule download_pgs_metadata:
    output:
        directory("data/raw/pgs_catalog/metadata/api")
    shell:
        """
        python scripts/download_pgs.py
        """


rule download_reactome:
    output:
        "data/raw/reactome/Ensembl2Reactome.txt",
        "data/raw/reactome/NCBI2Reactome.txt",
        "data/raw/reactome/Ensembl2Reactome_All_Levels.txt",
        "data/raw/reactome/UniProt2Reactome.txt",
    shell:
        "python scripts/download_reactome.py"
