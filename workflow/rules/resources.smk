rule get_genome_fasta:
    output:
        resources.fasta,
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log",
    conda:
        "../envs/bismark.yaml"
    threads: 1
    shell:
        "wget -q {params.url} -O {output}.gz 2> {log} ;"
        "pigz -d -c {output}.gz > {output}"


rule get_control_fasta:
    output:
        resources.control_fasta,
    params:
        url=resources.control_fasta_url,
    retries: 3
    log:
        "logs/resources/get_control_fasta.log",
    conda:
        "../envs/bismark.yaml"
    threads: 1
    shell:
        "wget -q {params.url} -O {output} 2> {log}"


rule combine_fasta:
    input:
        genome=resources.fasta,
        control=resources.control_fasta,
    output:
        "resources/combined_genome.fa",
    log:
        "logs/resources/combine_fasta.log",
    threads: 1
    conda:
        "../envs/bismark.yaml"
    shell:
        "cat {input.genome} {input.control}  > {output} 2> {log}"


rule index_fasta:
    input:
        "resources/combined_genome.fa",
    output:
        "resources/combined_genome.fa.fai",
    log:
        "logs/resources/index_fasta.log",
    threads: 1
    conda:
        "../envs/bismark.yaml"
    shell:
        "samtools faidx {input} 2> {log}"
    

rule chrom_sizes:
    input:
        "resources/combined_genome.fa.fai",
    output:
        "resources/chrom_sizes.txt",
    log:
        "logs/resources/chrom_sizes.log",
    threads: 1
    conda:
        "../envs/bismark.yaml"
    shell:
        "cut -f1,2 {input} > {output} 2> {log}"


use rule get_genome_fasta as get_gtf with:
    output:
        resources.gtf,
    params:
        url=resources.gtf_url,
    log:
        "logs/resources/get_gtf.log",


rule bismark_genome_preparation:
    input:
        fasta="resources/combined_genome.fa",
    output:
        directory("resources/Bisulfite_Genome")
    log:
        "logs/resources/bismark_genome_preparation.log"
    threads: 40 # make sure to assign half of this to bismark
    resources:
        runtime=360,
        mem_mb=60000,
    conda:
        "../envs/bismark.yaml"
    script:
        "../scripts/bismark_genome_preparation.py"

# Annotate CpGs in the genome
# -----------------------------------------------------
rule find_cpgs:
    input:
        resources.fasta,
    output:
        "resources/cpg_sites.bed",
    log:
        "logs/resources/find_cpgs.log",
    threads: 10
    resources:
        runtime=120,
        mem_mb=10000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "python workflow/scripts/find_cpgs.py {input} {output} {log}"

# Create CpG probe BED file
# -----------------------------------------------------
rule create_cpg_probes:
    input:
        "resources/cpg_sites.bed",
    output:
        "resources/cpg_probes.bed",
    params:
        # Number of CpGs per probe
        n=config["boxplot"]["cpg_n"],
    log:
        "logs/resources/create_cpg_probes.log",
    threads: 10
    resources:
        runtime=120,
        mem_mb=10000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "python workflow/scripts/create_cpg_probes.py {input} {output} {params.n} {log}"
