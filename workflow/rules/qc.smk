# make QC report
# -----------------------------------------------------
rule fastqc:
    input:
        fastq="results/simulate_reads/{sample}.bwa.{read}.fastq.gz",
    output:
        html="results/fastqc/{sample}.bwa.{read}_fastqc.html",
        zip="results/fastqc/{sample}.bwa.{read}_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "results/fastqc/{sample}.bwa.{read}.log",
    threads: 4
    wrapper:
        "v6.0.0/bio/fastqc"


# run multiQC on tool output
# -----------------------------------------------------
rule multiqc:
    input:
        expand(
            "results/fastqc/{sample}.bwa.{read}_fastqc.{ext}",
            sample=SAMPLES,
            read=["read1", "read2"],
            ext=["html", "zip"],
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    params:
        extra="--verbose --dirs",
    log:
        "results/multiqc/multiqc.log",
    wrapper:
        "v6.0.0/bio/multiqc"


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
        cpgot=expand(
            "results/bismark/{sample}/CpG_OT_{sample}.deduplicated.txt.gz",
            sample=SAMPLES,
        ),
        cpgob=expand(
            "results/bismark/{sample}/CpG_OB_{sample}.deduplicated.txt.gz",
            sample=SAMPLES,
        ),
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
        "sed 's/^\s*//;s/\s\s*/\t/g;s/$/\t{wildcards.sample}/' "
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