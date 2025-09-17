

rule get_fasta:
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


rule unpack_fasta:
    input:
        f"{resources.fasta}.gz",
    output:
        resources.fasta
    log:
        "logs/resources/unpack_fasta.log",
    conda:
        "../envs/bismark.yaml"
    threads: 1
    shell:
        "pigz -df {input} > {output} 2> {log}"

rule index_fasta:
    input:
        resources.fasta,
    output:
        f"{resources.fasta}.fai"
    log:
        "logs/resources/index_fasta.log",
    threads: 1
    conda:
        "../envs/bismark.yaml"
    shell:
        "samtools faidx {input} 2> {log}"
    

rule chrom_sizes:
    input:
        f"{resources.fasta}.fai"
    output:
        "resources/chrom_sizes.txt",
    log:
        "logs/resources/chrom_sizes.log",
    threads: 1
    conda:
        "../envs/bismark.yaml"
    shell:
        "cut -f1,2 {input} > {output} 2> {log}"


use rule get_fasta as get_gtf with:
    output:
        resources.gtf,
    params:
        url=resources.gtf_url,
    log:
        "logs/resources/get_gtf.log",


rule bismark_genome_preparation:
    input:
        fasta=resources.fasta,
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
