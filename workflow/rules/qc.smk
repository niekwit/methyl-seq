# Make QC report of raw data
# -----------------------------------------------------
rule fastqc:
    input:
        fastq="reads/{sample}{read}.fastq.gz",
    output:
        html="results/fastqc/{sample}{read}_fastqc.html",
        zip="results/fastqc/{sample}{read}_fastqc.zip",
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "results/fastqc/{sample}{read}.log",
    threads: 4
    resources:
        mem_mb = 1024,
        runtime = 30,
    wrapper:
        "v7.6.0/bio/fastqc"


# Run multiQC on FastQC output
# -----------------------------------------------------
rule multiqc:
    input:
        expand(
            "results/fastqc/{sample}{read}_fastqc.zip",
            sample=SAMPLES,
            read=["_R1_001", "_R2_001"] if PAIRED_END else [""],
            ext=["html", "zip"],
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    params:
        extra="--verbose --dirs",
    log:
        "logs/multiqc/fastqc.log",
    wrapper:
        "v8.1.1/bio/multiqc"


# Create control CpG coverage files
# -----------------------------------------------------
### AIM
# Count methylated/unmethylated CpGs in control DNA (lambda and pUC19)
# Output structure:
# count chrom pos methylation_status condition
# Example:
# 10  phage_lambda  123  Z  WT_1
# 9   plasmid_puc19c  1234  z  WT_1
# 1   plasmid_puc19c  1234  Z  WT_1

# NOTE: if a methylation call has no count for a probe
# it will not be represented in the output file (due to uniq -c).
# This will be corrected in R by adding 0 counts for missing methylation calls.
rule cpg_coverage_control_dna:
    input:
        cpgot="results/bismark/{sample}/CpG_OT_{sample}.deduplicated.txt.gz",
        cpgob="results/bismark/{sample}/CpG_OB_{sample}.deduplicated.txt.gz",
    output:
        cov="results/bismark/{sample}/CpG_{sample}.coverage.txt",
    log:
        "logs/cpg_coverage/{sample}.log",
    threads: 2
    resources:
        runtime=30,
    conda:
        "../envs/bismark.yaml"
    shell:
        "zcat {input.cpgot} {input.cpgob} | "
        "grep -E 'phage_lambda|plasmid_puc19c' | "
        "awk -v OFS=\"\t\" '{{print $3, $4, $5}}' | "
        "sort -k1,1 -k2,2n | "
        "uniq -c | "
        r"sed 's/^\s*//;s/\s\s*/\t/g;s/$/\t{wildcards.sample}/' "
        " > {output.cov} 2> {log}"


# Methylation converion rate calculation
# -----------------------------------------------------
rule plot_methylation_conversion_rate:
    input:
        cov=expand("results/bismark/{sample}/CpG_{sample}.coverage.txt", sample=SAMPLES),
    output:
        csv="results/plots/methylation_conversion_rate.csv",
        pdf="results/plots/methylation_conversion_rate.pdf",
    log:
        "logs/methylation_conversion_rate.log",
    threads: 2
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_methylation_conversion_rate.R"


if PAIRED_END:
    # Sort non-duplicated BAM files
    # -----------------------------------------------------
    rule sort_deduplicated_bam:
        input:
            "results/bismark/{sample}/{sample}_R1_bismark_bt2_pe.bam",
        output:
            temp("results/bismark/{sample}/{sample}.deduplicated.sorted.bam"),
        log:
            "logs/sort_deduplicated_bam/{sample}.log",
        threads: 4
        resources:
            runtime=45,
        wrapper:
            "v9.0.0/bio/samtools/sort"
else:
    # Sort non-duplicated BAM files
    # -----------------------------------------------------
    rule sort_deduplicated_bam:
        input:
            "results/bismark/{sample}/{sample}_bismark_bt2.bam",
        output:
            temp("results/bismark/{sample}/{sample}.deduplicated.sorted.bam"),
        log:
            "logs/sort_deduplicated_bam/{sample}.log",
        threads: 4
        resources:
            runtime=45,
        wrapper:
            "v9.0.0/bio/samtools/sort"

    # Analyse sequencing depth with preseq
    # -----------------------------------------------------
    rule preseq_lc_extrap_bam:
        input:
            "results/bismark/{sample}/{sample}.deduplicated.sorted.bam",
        output:
            "results/preseq/{sample}.txt",
        params:
            "-v"   #optional parameters
        log:
            "logs/preseq/{sample}.log",
        wrapper:
            "v2.10.0/bio/preseq/lc_extrap"

# Plot preseq results
# -----------------------------------------------------
rule plot_preseq:
    input:
        preseq=expand("results/preseq/{sample}.txt", sample=SAMPLES),
    output:
        pdf=expand("results/plots/library_complexity/{sample}.pdf", sample=SAMPLES),
        summary="results/preseq/library_complexity_summary.txt",
    log:
        "logs/plot_preseq.log",
    threads: 2
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_preseq.R"