# Convert CpG methylation calls of OT to BED (merge replicates)
# -----------------------------------------------------
rule OT_methylation_calls_to_bed:
    input:
        cpgot=lambda wildcards: expand(
            "results/bismark/{sample}/CpG_OT_{sample}.deduplicated.txt.gz",
            sample=[s for s in SAMPLES if s.startswith(wildcards.condition)],
        ),
    output:
        bed=temp("results/bed/CpG_OT_{condition}.bed"),
    log:
        "logs/methylation_calls_to_bed/OT_{condition}.log",
    threads: 4
    resources:
        runtime=240,
    conda:
        "../envs/bismark.yaml"
    shell:
        "zcat {input.cpgot} | "
        "grep -vE 'phage_lambda|plasmid_puc19c' | "
        'awk -v OFS="\t" \'{{print $3, $4, $4+1, $1, $5, "+"}}\' > {output.bed}'


# Convert CpG methylation calls of OB to BED (merge replicates)
# -----------------------------------------------------
rule OB_methylation_calls_to_bed:
    input:
        cpgob=lambda wildcards: expand(
            "results/bismark/{sample}/CpG_OB_{sample}.deduplicated.txt.gz",
            sample=[s for s in SAMPLES if s.startswith(wildcards.condition)],
        ),
    output:
        bed=temp("results/bed/CpG_OB_{condition}.bed"),
    log:
        "logs/methylation_calls_to_bed/OB_{condition}.log",
    threads: 4
    resources:
        runtime=240,
    conda:
        "../envs/bismark.yaml"
    shell:
        "zcat {input.cpgob} | "
        "grep -vE 'phage_lambda|plasmid_puc19c' | "
        'awk -v OFS="\t" \'{{print $3, $4, $4+1, $1, $5, "-"}}\' > {output.bed}'


# Merge CpG methylation calls of OT and OB to a single BED file
# -----------------------------------------------------
rule merge_strand_methylation_calls_to_bed:
    input:
        ot="results/bed/CpG_OT_{condition}.bed",
        ob="results/bed/CpG_OB_{condition}.bed",
    output:
        bed=temp("results/bed/CpG_merged_{condition}.bed"),
    log:
        "logs/methylation_calls_to_bed/merge_{condition}.log",
    threads: 2
    resources:
        runtime=60,
    conda:
        "../envs/bismark.yaml"
    shell:
        "cat {input.ot} {input.ob} | sort -k1,1 -k2,2n > {output.bed}"
