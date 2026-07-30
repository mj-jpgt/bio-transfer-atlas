"""Restartable, private MAG binning, consensus, QC, and taxonomy rules."""


BINNING = config["binning"]
DATABASES = config["databases"]
CHECKM2_DATABASE = os.environ.get("CHECKM2DB", DATABASES["checkm2"])
GUNC_DATABASE = os.environ.get("GUNC_DB", DATABASES["gunc"])
GTDBTK_DATABASE = os.environ.get("GTDBTK_DATA_PATH", DATABASES["gtdbtk"])
ASSEMBLY_INDEX_SUFFIXES = (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")


rule mag_consensus:
    input:
        expand(
            f"{WORK}/binning/{{sample}}/{{mode}}/dastool_scaffolds2bin.tsv",
            sample=SAMPLE_IDS,
            mode=FILTER_MODES,
        )


rule mag_endpoint_inputs:
    input:
        expand(
            f"{WORK}/binning/{{sample}}/{{mode}}/contig_bins.tsv",
            sample=SAMPLE_IDS,
            mode=FILTER_MODES,
        ),
        expand(
            f"{WORK}/binning/{{sample}}/{{mode}}/bin_qc.tsv",
            sample=SAMPLE_IDS,
            mode=FILTER_MODES,
        )


rule mag_database_preflight:
    output:
        ok=touch(f"{WORK}/databases/mag/preflight.ok"),
        versions=f"{WORK}/databases/mag/tool_versions.txt",
    log:
        "logs/mag/database_preflight.log",
    params:
        checkm2=CHECKM2_DATABASE,
        gunc=GUNC_DATABASE,
        gtdbtk=GTDBTK_DATABASE,
        release=DATABASES["gtdb_release"],
        temporary=f"{WORK}/databases/mag/.preflight.partial",
    threads:
        2
    resources:
        mem_mb=8000,
    conda:
        "../../envs/qc.yaml"
    shell:
        """
        mkdir -p $(dirname {log:q}) $(dirname {output.ok:q})
        rm -rf {params.temporary:q}
        mkdir -p {params.temporary:q}/gunc
        exec > {log:q} 2>&1
        test -s {params.checkm2:q}
        gunc check --db_file {params.gunc:q} \
          --out_dir {params.temporary:q}/gunc
        GTDBTK_DATA_PATH={params.gtdbtk:q} gtdbtk check_install
        {{
          checkm2 --version
          gunc --version
          GTDBTK_DATA_PATH={params.gtdbtk:q} gtdbtk --version
          printf 'GTDB release: %s\n' {params.release:q}
        }} > {output.versions:q}.partial
        mv {output.versions:q}.partial {output.versions:q}
        rm -rf {params.temporary:q}
        """


rule build_assembly_coverage_index:
    input:
        assembly=f"{WORK}/assembly/{{sample}}/{{mode}}/final.contigs.fa",
    output:
        index=[
            f"{WORK}/coverage/{{sample}}/{{mode}}/assembly{suffix}"
            for suffix in ASSEMBLY_INDEX_SUFFIXES
        ],
    log:
        "logs/mag/coverage_index/{sample}.{mode}.log",
    params:
        prefix=lambda wildcards, output: str(output.index[0]).removesuffix(".1.bt2"),
    wildcard_constraints:
        mode="|".join(FILTER_MODES),
    threads:
        config["resources"]["coverage"]["threads"]
    resources:
        mem_mb=config["resources"]["coverage"]["mem_mb"],
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        mkdir -p $(dirname {log:q}) $(dirname {params.prefix:q})
        bowtie2-build --threads {threads} {input.assembly:q} {params.prefix:q} \
          > {log:q} 2>&1
        """


rule map_reads_for_coverage:
    input:
        r1=f"{WORK}/filtered/{{sample}}/{{mode}}_R1.fastq.gz",
        r2=f"{WORK}/filtered/{{sample}}/{{mode}}_R2.fastq.gz",
        index=[
            f"{WORK}/coverage/{{sample}}/{{mode}}/assembly{suffix}"
            for suffix in ASSEMBLY_INDEX_SUFFIXES
        ],
    output:
        bam=f"{WORK}/coverage/{{sample}}/{{mode}}/reads.sorted.bam",
        bai=f"{WORK}/coverage/{{sample}}/{{mode}}/reads.sorted.bam.bai",
    log:
        "logs/mag/coverage/{sample}.{mode}.log",
    params:
        prefix=lambda wildcards, input: str(input.index[0]).removesuffix(".1.bt2"),
    wildcard_constraints:
        mode="|".join(FILTER_MODES),
    threads:
        config["resources"]["coverage"]["threads"]
    resources:
        mem_mb=config["resources"]["coverage"]["mem_mb"],
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        mkdir -p $(dirname {log:q})
        set -o pipefail
        bowtie2 --very-sensitive --no-unal -x {params.prefix:q} \
          -1 {input.r1:q} -2 {input.r2:q} -p {threads} 2> {log:q} \
          | samtools view -b -F 4 - \
          | samtools sort -@ {threads} -o {output.bam:q}.partial -
        mv {output.bam:q}.partial {output.bam:q}
        samtools index -@ {threads} {output.bam:q} {output.bai:q}.partial
        mv {output.bai:q}.partial {output.bai:q}
        """


rule summarize_contig_depth:
    input:
        bam=f"{WORK}/coverage/{{sample}}/{{mode}}/reads.sorted.bam",
        bai=f"{WORK}/coverage/{{sample}}/{{mode}}/reads.sorted.bam.bai",
    output:
        depth=f"{WORK}/coverage/{{sample}}/{{mode}}/contig_depth.tsv",
    log:
        "logs/mag/depth/{sample}.{mode}.log",
    threads:
        config["resources"]["coverage"]["threads"]
    resources:
        mem_mb=config["resources"]["coverage"]["mem_mb"],
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        mkdir -p $(dirname {log:q})
        jgi_summarize_bam_contig_depths --outputDepth {output.depth:q}.partial \
          {input.bam:q} > {log:q} 2>&1
        mv {output.depth:q}.partial {output.depth:q}
        """


rule maxbin_abundance:
    input:
        depth=f"{WORK}/coverage/{{sample}}/{{mode}}/contig_depth.tsv",
    output:
        abundance=f"{WORK}/coverage/{{sample}}/{{mode}}/maxbin_abundance.tsv",
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        PYTHONPATH=src python workflow/scripts/mag_translate.py abundance \
          --depth {input.depth:q} --output {output.abundance:q}
        """


rule metabat2_bins:
    input:
        assembly=f"{WORK}/assembly/{{sample}}/{{mode}}/final.contigs.fa",
        depth=f"{WORK}/coverage/{{sample}}/{{mode}}/contig_depth.tsv",
    output:
        bins=directory(f"{WORK}/binning/{{sample}}/{{mode}}/metabat2_bins"),
        mapping=f"{WORK}/binning/{{sample}}/{{mode}}/metabat2.scaffolds2bin.tsv",
    log:
        "logs/mag/metabat2/{sample}.{mode}.log",
    params:
        minimum=BINNING["minimum_contig_length"],
        temporary=lambda wildcards: (
            f"{WORK}/binning/{wildcards.sample}/{wildcards.mode}/.metabat2.partial"
        ),
    threads:
        config["resources"]["binner"]["threads"]
    resources:
        mem_mb=config["resources"]["binner"]["mem_mb"],
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        mkdir -p $(dirname {log:q})
        rm -rf {params.temporary:q}
        mkdir -p {params.temporary:q}
        metabat2 -i {input.assembly:q} -a {input.depth:q} \
          -o {params.temporary:q}/bin -m {params.minimum} -t {threads} \
          > {log:q} 2>&1
        PYTHONPATH=src python workflow/scripts/mag_translate.py bins-to-map \
          --bin-dir {params.temporary:q} \
          --bin-prefix metabat2. --output {output.mapping:q}
        mv {params.temporary:q} {output.bins:q}
        """


rule maxbin2_bins:
    input:
        assembly=f"{WORK}/assembly/{{sample}}/{{mode}}/final.contigs.fa",
        abundance=f"{WORK}/coverage/{{sample}}/{{mode}}/maxbin_abundance.tsv",
    output:
        bins=directory(f"{WORK}/binning/{{sample}}/{{mode}}/maxbin2_bins"),
        mapping=f"{WORK}/binning/{{sample}}/{{mode}}/maxbin2.scaffolds2bin.tsv",
    log:
        "logs/mag/maxbin2/{sample}.{mode}.log",
    params:
        minimum=BINNING["minimum_contig_length"],
        markerset=BINNING["maxbin_markerset"],
        temporary=lambda wildcards: (
            f"{WORK}/binning/{wildcards.sample}/{wildcards.mode}/.maxbin2.partial"
        ),
    threads:
        config["resources"]["binner"]["threads"]
    resources:
        mem_mb=config["resources"]["binner"]["mem_mb"],
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        mkdir -p $(dirname {log:q})
        rm -rf {params.temporary:q}
        mkdir -p {params.temporary:q}
        run_MaxBin.pl -contig {input.assembly:q} -abund {input.abundance:q} \
          -out {params.temporary:q}/bin -thread {threads} \
          -min_contig_length {params.minimum} -markerset {params.markerset} \
          > {log:q} 2>&1
        PYTHONPATH=src python workflow/scripts/mag_translate.py bins-to-map \
          --bin-dir {params.temporary:q} \
          --bin-prefix maxbin2. --output {output.mapping:q}
        mv {params.temporary:q} {output.bins:q}
        """


rule concoct_bins:
    input:
        assembly=f"{WORK}/assembly/{{sample}}/{{mode}}/final.contigs.fa",
        bam=f"{WORK}/coverage/{{sample}}/{{mode}}/reads.sorted.bam",
        bai=f"{WORK}/coverage/{{sample}}/{{mode}}/reads.sorted.bam.bai",
    output:
        bins=directory(f"{WORK}/binning/{{sample}}/{{mode}}/concoct_bins"),
        mapping=f"{WORK}/binning/{{sample}}/{{mode}}/concoct.scaffolds2bin.tsv",
    log:
        "logs/mag/concoct/{sample}.{mode}.log",
    params:
        chunk=BINNING["concoct_chunk_length"],
        temporary=lambda wildcards: (
            f"{WORK}/binning/{wildcards.sample}/{wildcards.mode}/.concoct.partial"
        ),
    threads:
        config["resources"]["binner"]["threads"]
    resources:
        mem_mb=config["resources"]["binner"]["mem_mb"],
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        mkdir -p $(dirname {log:q})
        rm -rf {params.temporary:q}
        mkdir -p {params.temporary:q}/run {params.temporary:q}/bins
        exec > {log:q} 2>&1
        cut_up_fasta.py {input.assembly:q} -c {params.chunk} -o 0 --merge_last \
          -b {params.temporary:q}/cut.bed > {params.temporary:q}/cut.fa
        concoct_coverage_table.py {params.temporary:q}/cut.bed {input.bam:q} \
          > {params.temporary:q}/coverage.tsv
        concoct --threads {threads} \
          --composition_file {params.temporary:q}/cut.fa \
          --coverage_file {params.temporary:q}/coverage.tsv \
          -b {params.temporary:q}/run/
        merge_cutup_clustering.py \
          {params.temporary:q}/run/clustering_gt1000.csv \
          > {params.temporary:q}/clustering_merged.csv
        extract_fasta_bins.py {input.assembly:q} \
          {params.temporary:q}/clustering_merged.csv \
          --output_path {params.temporary:q}/bins
        PYTHONPATH=src python workflow/scripts/mag_translate.py bins-to-map \
          --bin-dir {params.temporary:q}/bins \
          --bin-prefix concoct. --output {output.mapping:q}
        mv {params.temporary:q}/bins {output.bins:q}
        rm -rf {params.temporary:q}
        """


rule dastool_consensus:
    input:
        assembly=f"{WORK}/assembly/{{sample}}/{{mode}}/final.contigs.fa",
        metabat=f"{WORK}/binning/{{sample}}/{{mode}}/metabat2.scaffolds2bin.tsv",
        maxbin=f"{WORK}/binning/{{sample}}/{{mode}}/maxbin2.scaffolds2bin.tsv",
        concoct=f"{WORK}/binning/{{sample}}/{{mode}}/concoct.scaffolds2bin.tsv",
    output:
        bins=directory(f"{WORK}/binning/{{sample}}/{{mode}}/dastool_bins"),
        mapping=f"{WORK}/binning/{{sample}}/{{mode}}/dastool_scaffolds2bin.tsv",
        summary=f"{WORK}/binning/{{sample}}/{{mode}}/dastool_summary.tsv",
    log:
        "logs/mag/dastool/{sample}.{mode}.log",
    params:
        score=BINNING["dastool_score_threshold"],
        engine=BINNING["dastool_search_engine"],
        temporary=lambda wildcards: (
            f"{WORK}/binning/{wildcards.sample}/{wildcards.mode}/.dastool.partial"
        ),
    threads:
        config["resources"]["dastool"]["threads"]
    resources:
        mem_mb=config["resources"]["dastool"]["mem_mb"],
    conda:
        "../../envs/assembly.yaml"
    shell:
        """
        mkdir -p $(dirname {log:q})
        rm -rf {params.temporary:q}
        mkdir -p {params.temporary:q}
        DAS_Tool \
          -i {input.metabat:q},{input.maxbin:q},{input.concoct:q} \
          -l metabat2,maxbin2,concoct -c {input.assembly:q} \
          -o {params.temporary:q}/consensus --threads {threads} \
          --search_engine {params.engine:q} --score_threshold {params.score} \
          --write_bins 1 --write_bin_evals 1 --create_plots 0 \
          > {log:q} 2>&1
        mv {params.temporary:q}/consensus_DASTool_bins {output.bins:q}
        mv {params.temporary:q}/consensus_DASTool_scaffolds2bin.txt \
          {output.mapping:q}
        mv {params.temporary:q}/consensus_DASTool_summary.txt {output.summary:q}
        rm -rf {params.temporary:q}
        """


rule checkm2_qc:
    input:
        bins=f"{WORK}/binning/{{sample}}/{{mode}}/dastool_bins",
        preflight=f"{WORK}/databases/mag/preflight.ok",
    output:
        report=f"{WORK}/binning/{{sample}}/{{mode}}/qc/checkm2_quality.tsv",
    log:
        "logs/mag/checkm2/{sample}.{mode}.log",
    params:
        database=CHECKM2_DATABASE,
        temporary=lambda wildcards: (
            f"{WORK}/binning/{wildcards.sample}/{wildcards.mode}/qc/"
            ".checkm2.partial"
        ),
    threads:
        config["resources"]["checkm2"]["threads"]
    resources:
        mem_mb=config["resources"]["checkm2"]["mem_mb"],
    conda:
        "../../envs/qc.yaml"
    shell:
        """
        mkdir -p $(dirname {log:q}) $(dirname {output.report:q})
        rm -rf {params.temporary:q}
        checkm2 predict --threads {threads} --input {input.bins:q} \
          --output-directory {params.temporary:q} \
          --database_path {params.database:q} > {log:q} 2>&1
        mv {params.temporary:q}/quality_report.tsv {output.report:q}
        rm -rf {params.temporary:q}
        """


rule gunc_qc:
    input:
        bins=f"{WORK}/binning/{{sample}}/{{mode}}/dastool_bins",
        preflight=f"{WORK}/databases/mag/preflight.ok",
    output:
        report=f"{WORK}/binning/{{sample}}/{{mode}}/qc/gunc_maxcss.tsv",
    log:
        "logs/mag/gunc/{sample}.{mode}.log",
    params:
        database=GUNC_DATABASE,
        temporary=lambda wildcards: (
            f"{WORK}/binning/{wildcards.sample}/{wildcards.mode}/qc/.gunc.partial"
        ),
    threads:
        config["resources"]["gunc"]["threads"]
    resources:
        mem_mb=config["resources"]["gunc"]["mem_mb"],
    conda:
        "../../envs/qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log:q}) $(dirname {output.report:q})
        rm -rf {params.temporary:q}
        mkdir -p {params.temporary:q}
        gunc run --input_dir {input.bins:q} --file_suffix .fa \
          --db_file {params.database:q} --threads {threads} \
          --out_dir {params.temporary:q} > {log:q} 2>&1
        reports=({params.temporary:q}/GUNC.*.maxCSS_level.tsv)
        test "${{#reports[@]}}" -eq 1
        mv "${{reports[0]}}" {output.report:q}
        rm -rf {params.temporary:q}
        """


rule gtdbtk_r220_taxonomy:
    input:
        bins=f"{WORK}/binning/{{sample}}/{{mode}}/dastool_bins",
        preflight=f"{WORK}/databases/mag/preflight.ok",
    output:
        bacterial=f"{WORK}/binning/{{sample}}/{{mode}}/taxonomy/gtdbtk.bac120.summary.tsv",
        archaeal=f"{WORK}/binning/{{sample}}/{{mode}}/taxonomy/gtdbtk.ar53.summary.tsv",
    log:
        "logs/mag/gtdbtk/{sample}.{mode}.log",
    params:
        database=GTDBTK_DATABASE,
        temporary=lambda wildcards: (
            f"{WORK}/binning/{wildcards.sample}/{wildcards.mode}/taxonomy/"
            ".gtdbtk.partial"
        ),
    threads:
        config["resources"]["gtdbtk"]["threads"]
    resources:
        mem_mb=config["resources"]["gtdbtk"]["mem_mb"],
    conda:
        "../../envs/qc.yaml"
    shell:
        """
        mkdir -p $(dirname {log:q}) $(dirname {output.bacterial:q})
        rm -rf {params.temporary:q}
        GTDBTK_DATA_PATH={params.database:q} gtdbtk classify_wf \
          --genome_dir {input.bins:q} --out_dir {params.temporary:q} \
          --extension fa --cpus {threads} --prefix gtdbtk \
          > {log:q} 2>&1
        mv {params.temporary:q}/gtdbtk.bac120.summary.tsv \
          {output.bacterial:q}
        mv {params.temporary:q}/gtdbtk.ar53.summary.tsv \
          {output.archaeal:q}
        rm -rf {params.temporary:q}
        """


rule mag_private_contract:
    input:
        dastool=f"{WORK}/binning/{{sample}}/{{mode}}/dastool_scaffolds2bin.tsv",
        checkm2=f"{WORK}/binning/{{sample}}/{{mode}}/qc/checkm2_quality.tsv",
        gunc=f"{WORK}/binning/{{sample}}/{{mode}}/qc/gunc_maxcss.tsv",
        bacterial=f"{WORK}/binning/{{sample}}/{{mode}}/taxonomy/gtdbtk.bac120.summary.tsv",
        archaeal=f"{WORK}/binning/{{sample}}/{{mode}}/taxonomy/gtdbtk.ar53.summary.tsv",
    output:
        bins=f"{WORK}/binning/{{sample}}/{{mode}}/contig_bins.tsv",
        qc=f"{WORK}/binning/{{sample}}/{{mode}}/bin_qc.tsv",
    log:
        "logs/mag/contract/{sample}.{mode}.log",
    threads:
        config["resources"]["mag_contract"]["threads"]
    resources:
        mem_mb=config["resources"]["mag_contract"]["mem_mb"],
    conda:
        "../../envs/qc.yaml"
    shell:
        """
        mkdir -p $(dirname {log:q})
        PYTHONPATH=src python workflow/scripts/mag_translate.py contract \
          --sample-id {wildcards.sample:q} --dastool-map {input.dastool:q} \
          --checkm2-report {input.checkm2:q} --gunc-report {input.gunc:q} \
          --gtdb-bacterial-summary {input.bacterial:q} \
          --gtdb-archaeal-summary {input.archaeal:q} \
          --contig-bins-output {output.bins:q} \
          --bin-qc-output {output.qc:q} > {log:q} 2>&1
        """
