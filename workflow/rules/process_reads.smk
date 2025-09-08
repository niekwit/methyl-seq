# Trim reads
# -----------------------------------------------------
rule trim_galore:
    input:
        ["reads/{sample}_R1_001.fastq.gz", "reads/{sample}_R2_001.fastq.gz"],
    output:
        fasta_fwd=temp("results/trimmed/{sample}_R1.fq.gz"),
        report_fwd="logs/trim_galore/{sample}_R1_trimming_report.txt",
        fasta_rev=temp("results/trimmed/{sample}_R2.fq.gz"),
        report_rev="logs/trim_galore/{sample}_R2_trimming_report.txt",
    threads: 4
    resources: 
        runtime=30,
    params:
        extra="--illumina -q 20",
    log:
        "logs/trim_galore/{sample}.log",
    wrapper:
        "v7.2.0/bio/trim_galore/pe"

# Align reads with Bismark
# -----------------------------------------------------
rule align:
    input:
        dir="resources/Bisulfite_Genome",
        r1="results/trimmed/{sample}_R1.fq.gz",
        r2="results/trimmed/{sample}_R2.fq.gz",
    output:
        bam=temp("results/bismark/{sample}_R1_bismark_bt2_pe.bam"),
    log:
        "logs/bismark_align/{sample}.log",
    threads: 12
    resources:
        runtime=60,
    conda:
        "../envs/bismark.yaml"
    shell:
        "bismark "
        "--genome resources/ "
        "-p {threads} "
        "-1 {input.r1} "
        "-2 {input.r2} "
        "-o results/aligned/ "
        "2> {log}"

# Deduplicate aligned reads with Bismark
# -----------------------------------------------------
rule deduplication:
    input:
        bam="results/bismark/{sample}_R1_bismark_bt2_pe.bam",
    output:
        bam="results/bismark/{sample}.deduplicated.bam",
    log:
        "logs/deduplication/{sample}.log",
    threads: 4
    resources:
        runtime=30,
    conda:
        "../envs/bismark.yaml"
    shell:
        "deduplicate_bismark "
        "--paired "
        "--outfile {wildcards.sample} "
        "--output_dir results/deduplicated/ "
        "--bam "
        "{input.bam} "
        "2> {log}"

# Extract methylation call for every single C analysed with Bismark
# -----------------------------------------------------
'''
rule methylation_extraction:
    input:
        bam="results/bismark/{sample}.deduplicated.bam",
    output:
        dir=directory("results/bismark/{sample}/")
    log:
        "logs/methylation_extraction/{sample}.log"
    threads: 4
    resources:
        runtime=30,
    conda:
        "../envs/bismark.yaml"
    shell:
        "bismark_methylation_extractor "
        "--paired-end "
        "--output_dir {output.dir} "
        "--bam "
        "--cytosine_report "
        "--gzip "
        "--no_header "
        "--buffer_size 10G "
        "--genome_folder resources/ "
        "{input.bam} "
        "2> {log}"
'''
# Extract nucleotide coverage
# -----------------------------------------------------
rule nucleotide_coverage:
    input:
        bam="results/bismark/{sample}.deduplicated.bam",
    output:
        stats="results/bismark/{sample}.deduplicated.nucleotide_stats.txt",
    params:
        dir=lambda wc, output: os.path.dirname(output.stats)
    log:
        "logs/nucleotide_coverage/{sample}.log"
    threads: 4
    resources:
        runtime=30,
    conda:
        "../envs/bismark.yaml"
    shell:
        "mkdir -p {params.dir}; "
        "bam2nuc "
        "--dir {params.dir} "
        "--genome_folder resources/ "
        "{input.bam} "
        "2> {log}"

'''
rule summary_report:
    input:
        expand("results/bismark/{sample}.bam", sample=SAMPLES),
    output:
        "results/bismark/report.html",
    log:
        "logs/bismark/summary_report.log"
    threads: 2
    resources:
        runtime=30,
    conda:
        "../envs/bismark.yaml"
    shell:
        "bismark2summary -o {output} {input}"
'''