# Filter CpG probes to keep only those covered by >10 reads (per condition)
# -----------------------------------------------------
rule filter_cpg_probes_for_reads:
    input:
        cpg="resources/cpg_probes.bed",
        meth="results/bed/CpG_merged_{condition}.bed",
    output:
        temp("results/bed/CpG_probes_{condition}_filtered.txt"),
    log:
        "logs/score_methylation_calls/{condition}.log",
    threads: 4
    resources:
        runtime=180,
        mem_mb=6000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bedtools intersect -wa -wb -sorted -a {input.cpg} -b {input.meth} | "
        "cut -f4,8 | "
        "sort -k 1,1 -k2,2n | "
        "uniq | "
        "cut -f1 | "
        "uniq -c | "
        "awk '$1 > 10 {{print $2}}' "
        "> {output}"


# Only keep CpG probes that are covered by >10 reads in ALL conditions
# -----------------------------------------------------
rule filter_cpg_probes_all_conditions:
    input:
        probes="resources/cpg_probes.bed",
        fprobes=expand(
            "results/bed/CpG_probes_{condition}_filtered.txt", condition=CONDITIONS
        ),
    output:
        "resources/filtered_cpg_probes.bed",
    log:
        "logs/score_methylation_calls/filter_all_conditions.log",
    threads: 4
    resources:
        runtime=120,
        mem_mb=6000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "sort {input.fprobes} | "
        "uniq -d | "
        "grep -wF -f - {input.probes} | "
        "sort -k1,1 -k2,2n "
        "> {output}"


# Filter CpG probes to keep only those in specified regions
# -----------------------------------------------------
rule filter_cpg_probes_for_regions:
    input:
        probes="resources/filtered_cpg_probes.bed",
        regions="bed/{region}.bed",
    output:
        region_probes=temp("resources/cpg_probes_{region}.bed"),
    log:
        "logs/resources/filter_probes_{region}.log",
    threads: 1
    resources:
        runtime=15,
        mem_mb=2000,
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/filter_cpg_probes_for_regions.sh"


# Create input for boxplot
# -----------------------------------------------------
rule boxplot_data:
    input:
        probes="resources/cpg_probes_{region}.bed",
        meth="results/bed/CpG_merged_{condition}.bed",
    output:
        temp("results/boxplot/CpG_methylation_{condition}_{region}.txt"),
    log:
        "logs/score_methylation_calls/boxplot_{condition}_{region}.log",
    threads: 1
    resources:
        runtime=360,
        mem_mb=2000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bedtools intersect -sorted -wa -wb -a {input.probes} -b {input.meth} | "
        "cut -f4,9 | "
        "sort -k1,1 -k2,2n | "
        "uniq -c | "
        "sed 's/$/\t{wildcards.region}\t{wildcards.condition}/' > {output}"


# Combine boxplot data for all conditions and regions
# -----------------------------------------------------
rule combine_boxplot_data:
    input:
        data=expand(
            "results/boxplot/CpG_methylation_{condition}_{region}.txt",
            condition=CONDITIONS,
            region=REGIONS,
        ),
    output:
        "results/boxplot/CpG_methylation_all_conditions_all_regions.txt",
    log:
        "logs/boxplot/combine_data.log",
    threads: 1
    resources:
        runtime=10,
        mem_mb=2000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "cat {input.data} | "
        r"sed 's/^\s*//;s/\s/\t/g' > {output}"


# Plot CpG methylation boxplot
# -----------------------------------------------------
rule plot_boxplot:
    input:
        "results/boxplot/CpG_methylation_all_conditions_all_regions.txt",
    output:
        pdf="results/plots/boxplots.pdf",
    params:
        regions=REGIONS,
    log:
        "logs/boxplot/plot_boxplots.log",
    threads: 1
    resources:
        runtime=30,
        mem_mb=2000,
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_boxplot.R"
