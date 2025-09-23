
if PAIRED_END:
    # Trim reads
    # -----------------------------------------------------
    rule trim_galore_pe:
        input:
            ["reads/{sample}_R1_001.fastq.gz", "reads/{sample}_R2_001.fastq.gz"],
        output:
            fasta_fwd=temp("results/trimmed/{sample}_R1.fq.gz"),
            report_fwd="logs/trim_galore/{sample}_R1_trimming_report.txt",
            fasta_rev=temp("results/trimmed/{sample}_R2.fq.gz"),
            report_rev="logs/trim_galore/{sample}_R2_trimming_report.txt",
        threads: 4
        resources: 
            runtime=120,
            mem_mb=20000,
            tmpdir=config["temp_dir"],
        params:
            extra=f"--illumina -q 20 {config["trim_galore_args"]}",
        log:
            "logs/trim_galore_pe/{sample}.log",
        wrapper:
            "v7.2.0/bio/trim_galore/pe"

    # Align reads with Bismark
    # -----------------------------------------------------
    rule align_pe:
        input:
            dir="resources/Bisulfite_Genome",
            r1="results/trimmed/{sample}_R1.fq.gz",
            r2="results/trimmed/{sample}_R2.fq.gz",
        output:
            bam=temp("results/bismark/{sample}/{sample}_R1_bismark_bt2_pe.bam"),
        params:
            outdir=lambda wc, output: os.path.dirname(output.bam),
        log:
            "logs/bismark_align/{sample}.log",
        threads: 12
        resources:
            runtime=2000,
        conda:
            "../envs/bismark.yaml"
        shell:
            "bismark "
            "--genome resources/ "
            "-p {threads} "
            "-1 {input.r1} "
            "-2 {input.r2} "
            "-o {params.outdir} "
            "2> {log}"
else:
    # Trim reads
    # -----------------------------------------------------
    rule trim_galore_se:
        input:
            "reads/{sample}.fastq.gz",
        output:
            fasta=temp("results/trimmed/{sample}.fq.gz"),
            report=temp("results/trimmed/{sample}_report.txt"),
        params:
            extra=f"--illumina -q 20 {config["trim_galore_args"]}",
        log:
            "logs/trim_galore_se/{sample}.log",
        wrapper:
            "v7.2.0/bio/trim_galore/se"

    # Align reads with Bismark
    # -----------------------------------------------------
    rule align_se:
        input:
            dir="resources/Bisulfite_Genome",
            fq="results/trimmed/{sample}.fq.gz",
        output:
            bam=temp("results/bismark/{sample}/{sample}_bismark_bt2.bam"),
        params:
            outdir=lambda wc, output: os.path.dirname(output.bam),
        log:
            "logs/bismark_align/{sample}.log",
        threads: 12
        resources:
            runtime=2000,
        conda:
            "../envs/bismark.yaml"
        shell:
            "bismark "
            "-o {params.outdir} "
            "--genome resources/ "
            "-p {threads} "
            "{input.fq} "
            "2> {log}"

# Deduplicate aligned reads with Bismark
# -----------------------------------------------------
rule deduplication:
    input:
        unpack(dedup_input)
    output:
        bam="results/bismark/{sample}/{sample}.deduplicated.bam",
    params:
        outdir=lambda wc, output: os.path.dirname(output.bam),
    log:
        "logs/deduplication/{sample}.log",
    threads: 4
    resources:
        runtime=600,
    conda:
        "../envs/bismark.yaml"
    shell:
        "deduplicate_bismark "
        "--paired "
        "--outfile {wildcards.sample} "
        "--output_dir {params.outdir} "
        "--bam "
        "{input.bam} "
        "2> {log}"

# Extract methylation call for every single C analysed with Bismark
# -----------------------------------------------------
rule methylation_extraction:
    input:
        bam="results/bismark/{sample}/{sample}.deduplicated.bam",
    output:
        sreport="results/bismark/{sample}/{sample}.deduplicated_splitting_report.txt",
        mbias="results/bismark/{sample}/{sample}.deduplicated.M-bias.txt",
        bg="results/bismark/{sample}/{sample}.deduplicated.bedGraph.gz", #CpG only
        cpgot="results/bismark/{sample}/CpG_OT_{sample}.deduplicated.txt.gz",
        cpgob="results/bismark/{sample}/CpG_OB_{sample}.deduplicated.txt.gz",
    params:
        outdir=lambda wc, output: os.path.dirname(output.sreport),
    log:
        "logs/methylation_extraction/{sample}.log"
    threads: 4
    resources:
        runtime=600,
    conda:
        "../envs/bismark.yaml"
    shell:
        "bismark_methylation_extractor "
        "--paired-end "
        "--no_overlap "
        "--output_dir {params.outdir} "
        "--bedgraph "
        "--cytosine_report "
        "--gzip "
        "--no_header "
        "--buffer_size 10G "
        "--genome_folder resources/ "
        "--multicore {threads} "
        "{input.bam} "
        "2> {log}"

# Extract nucleotide coverage
# -----------------------------------------------------
rule nucleotide_coverage:
    input:
        bam="results/bismark/{sample}/{sample}.deduplicated.bam",
    output:
        stats="results/bismark/{sample}/{sample}.deduplicated.nucleotide_stats.txt",
    params:
        dir=lambda wc, output: os.path.dirname(output.stats)
    log:
        "logs/nucleotide_coverage/{sample}.log"
    threads: 4
    resources:
        runtime=600,
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
        expand("results/bismark/{sample}/{sample}.bam", sample=SAMPLES),
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