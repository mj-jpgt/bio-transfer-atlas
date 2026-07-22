"""
Snakemake rules: VCF -> PLINK2 pgen, QC, LD-pruning, PCA.
"""

VCF_TEMPLATE = "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"


rule vcf_to_pgen:
    input:
        vcf = "data/raw/1000g/vcf/" + VCF_TEMPLATE,
    output:
        multiext("data/interim/1000g/chr{chrom}", ".pgen", ".psam", ".pvar"),
    shell:
        """
        plink2 \
          --vcf {input.vcf} \
          --make-pgen \
          --out data/interim/1000g/chr{wildcards.chrom}
        """


rule genotype_qc:
    input:
        multiext("data/interim/1000g/chr{chrom}", ".pgen", ".psam", ".pvar"),
    output:
        multiext("data/interim/1000g/chr{chrom}.qc", ".pgen", ".psam", ".pvar"),
    shell:
        """
        plink2 \
          --pfile data/interim/1000g/chr{wildcards.chrom} \
          --geno 0.02 \
          --mind 0.02 \
          --maf 0.01 \
          --make-pgen \
          --out data/interim/1000g/chr{wildcards.chrom}.qc
        """


rule ld_prune:
    input:
        multiext("data/interim/1000g/chr{chrom}.qc", ".pgen", ".psam", ".pvar"),
    output:
        prune_in  = "data/interim/1000g/chr{chrom}.prune.in",
        prune_out = "data/interim/1000g/chr{chrom}.prune.out",
    shell:
        """
        plink2 \
          --pfile data/interim/1000g/chr{wildcards.chrom}.qc \
          --indep-pairwise 200 50 0.2 \
          --out data/interim/1000g/chr{wildcards.chrom}.prune
        """


rule compute_pca:
    input:
        pgen      = "data/interim/1000g/chr{chrom}.qc.pgen",
        prune_in  = "data/interim/1000g/chr{chrom}.prune.in",
    output:
        eigenvec  = "data/processed/1000g/chr{chrom}_pca.eigenvec",
        eigenval  = "data/processed/1000g/chr{chrom}_pca.eigenval",
    shell:
        """
        mkdir -p data/processed/1000g
        plink2 \
          --pfile data/interim/1000g/chr{wildcards.chrom}.qc \
          --extract {input.prune_in} \
          --pca 20 \
          --out data/processed/1000g/chr{wildcards.chrom}_pca
        """
