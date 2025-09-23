rule get_genome_fasta:
    output:
        f"{resources.fasta}.gz",
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log",
    conda:
        "../envs/bismark.yaml"
    threads: 1
    shell:
        "wget -q {params.url} -O {output} 2> {log}"


use rule get_genome_fasta as get_control_fasta with:
    output:
        f"{resources.control_fasta}",
    params:
        url=resources.control_fasta_url,
    log:
        "logs/resources/get_control_fasta.log",


rule combine_fasta:
    input:
        genome=f"{resources.fasta}.gz",
        control=resources.control_fasta,
    output:
        "resources/combined_genome.fa",
    log:
        "logs/resources/combine_fasta.log",
    threads: 1
    conda:
        "../envs/bismark.yaml"
    shell:
        "cat <(zcat {input.genome}) {input.control}  > {output} 2> {log}"


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
    conda:
        "../envs/bismark.yaml"
    script:
        "../scripts/bismark_genome_preparation.py"
