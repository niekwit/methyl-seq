if PAIRED_END:

    # Trim reads
    # -----------------------------------------------------
    rule trim_galore_pe:
        input:
            ["reads/{sample}_R1_001.fastq.gz", "reads/{sample}_R2_001.fastq.gz"],
        output:
            fasta_fwd=temp("results/trimmed/{sample}_R1.fq.gz"),
            report_fwd="logs/trim_galore_pe/{sample}_R1_trimming_report.txt",
            fasta_rev=temp("results/trimmed/{sample}_R2.fq.gz"),
            report_rev="logs/trim_galore_pe/{sample}_R2_trimming_report.txt",
        threads: 4
        resources:
            runtime=180,
            mem_mb=20000,
            tmpdir=config["temp_dir"],
        params:
            extra=f"--illumina -q 20 {config["trim_galore_args"]}",
        log:
            "logs/trim_galore_pe/{sample}.log",
        conda:
            "../envs/trim_galore.yml"
        wrapper:
            "v9.9.0/bio/trim_galore/pe"

    # Align reads with Bismark
    # -----------------------------------------------------
    rule align_pe:
        input:
            dir="resources/Bisulfite_Genome",
            # Not referenced directly in the shell command below: bismark's
            # aligner does its own <*.fa> scan of --genome (resources/) at
            # runtime to load genomic sequence for methylation calling, a
            # dependency Snakemake can't see from "dir" alone. Without this,
            # --all-temp (used in CI) can delete combined_genome.fa right
            # after bismark_genome_preparation finishes with it -- before
            # this rule runs -- and the aligner then silently falls back to
            # whatever other *.fa/*.fa.gz happens to sit in resources/ (e.g.
            # the raw, uncombined genome.fa.gz), loading only that genome
            # and discarding every read outside it ("genomic sequence could
            # not be extracted").
            fasta="resources/combined_genome.fa",
            r1="results/trimmed/{sample}_R1.fq.gz",
            r2="results/trimmed/{sample}_R2.fq.gz",
        output:
            bam="results/bismark/{sample}/{sample}_R1_bismark_bt2_pe.bam",
        params:
            outdir=lambda wc, output: os.path.dirname(output.bam),
            extra=config["bismark"]["align"],
        log:
            "logs/bismark_align/{sample}.log",
        threads: 12
        resources:
            runtime=2000,
            mem_mb=60000,
        conda:
            "../envs/bismark.yaml"
        shell:
            "bismark "
            "--genome resources/ "
            "-p {threads} "
            "-1 {input.r1} "
            "-2 {input.r2} "
            "-o {params.outdir} "
            "{params.extra} "
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
        resources:
            runtime=180,
            mem_mb=20000,
            tmpdir=config["temp_dir"],
        log:
            "logs/trim_galore_se/{sample}.log",
        conda:
            "../envs/trim_galore.yml"
        wrapper:
            "v9.9.0/bio/trim_galore/se"

    # Align reads with Bismark
    # -----------------------------------------------------
    rule align_se:
        input:
            dir="resources/Bisulfite_Genome",
            # See the matching comment in align_pe above.
            fasta="resources/combined_genome.fa",
            fq="results/trimmed/{sample}.fq.gz",
        output:
            bam="results/bismark/{sample}/{sample}_bismark_bt2.bam",
        params:
            outdir=lambda wc, output: os.path.dirname(output.bam),
            extra=config["bismark"]["align"],
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
            "{params.extra} "
            "2> {log}"


# Deduplicate aligned reads with Bismark
# -----------------------------------------------------
rule deduplication:
    input:
        unpack(dedup_input),
    output:
        bam="results/bismark/{sample}/{sample}.deduplicated.bam",
    params:
        outdir=lambda wc, output: os.path.dirname(output.bam),
        paired="--paired" if PAIRED_END else "",
        extra=config["bismark"]["deduplicate"],
    log:
        "logs/deduplication/{sample}.log",
    threads: 4
    resources:
        runtime=600,
    conda:
        "../envs/bismark.yaml"
    shell:
        "deduplicate_bismark "
        "{params.paired} "
        "--outfile {wildcards.sample} "
        "--output_dir {params.outdir} "
        "--bam "
        "{input.bam} "
        "{params.extra} "
        "2> {log}"


# Extract methylation call for every single C analysed with Bismark
# -----------------------------------------------------
rule methylation_extraction:
    input:
        bam="results/bismark/{sample}/{sample}.deduplicated.bam",
        # See the matching comment in align_pe (process_reads.smk):
        # --cytosine_report reads --genome_folder's *.fa directly too, a
        # dependency Snakemake can't otherwise see from "bam" alone.
        fasta="resources/combined_genome.fa",
    output:
        sreport="results/bismark/{sample}/{sample}.deduplicated_splitting_report.txt",
        mbias="results/bismark/{sample}/{sample}.deduplicated.M-bias.txt",
        bg="results/bismark/{sample}/{sample}.deduplicated.bedGraph.gz",  #CpG only
        cpgot="results/bismark/{sample}/CpG_OT_{sample}.deduplicated.txt.gz",
        cpgob="results/bismark/{sample}/CpG_OB_{sample}.deduplicated.txt.gz",
        cov="results/bismark/{sample}/{sample}.deduplicated.bismark.cov.gz",
    params:
        outdir=lambda wc, output: os.path.dirname(output.sreport),
        genome_abs=os.path.abspath("resources/"),
        paired="--paired-end" if PAIRED_END else "",
        # bismark_methylation_extractor handles both extraction and (via
        # --cytosine_report) coverage/cytosine-report generation in one
        # invocation -- there's no separate coverage2cytosine step in this
        # workflow -- so config bismark:extract and bismark:coverage both
        # apply here.
        extra_extract=config["bismark"]["extract"],
        extra_coverage=config["bismark"]["coverage"],
    log:
        "logs/methylation_extraction/{sample}.log",
    threads: 4
    resources:
        runtime=600,
    conda:
        "../envs/bismark.yaml"
    shell:
        "bismark_methylation_extractor "
        "{params.paired} "
        "--no_overlap "
        "--output_dir {params.outdir} "
        "--bedgraph "
        "--cytosine_report "
        "--gzip "
        "--no_header "
        "--buffer_size 10G "
        "--genome_folder {params.genome_abs} "
        "--multicore {threads} "
        "{input.bam} "
        "{params.extra_extract} "
        "{params.extra_coverage} "
        "2> {log}"


# Extract nucleotide coverage
# -----------------------------------------------------
rule nucleotide_coverage:
    input:
        bam="results/bismark/{sample}/{sample}.deduplicated.bam",
        # See the matching comment in align_pe (process_reads.smk):
        # bam2nuc reads --genome_folder's *.fa directly too.
        fasta="resources/combined_genome.fa",
    output:
        stats="results/bismark/{sample}/{sample}.deduplicated.nucleotide_stats.txt",
    params:
        dir=lambda wc, output: os.path.dirname(output.stats),
    log:
        "logs/nucleotide_coverage/{sample}.log",
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


rule multiqc_bismark:
    input:
        expand(
            "results/bismark/{sample}/{sample}.deduplicated.nucleotide_stats.txt",
            sample=SAMPLES,
        ),
    output:
        "results/multiqc/multiqc_bismark.html",
    log:
        "logs/multiqc/bismark.log",
    threads: 2
    resources:
        runtime=30,
    wrapper:
        "v8.1.1/bio/multiqc"


"""
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
"""
