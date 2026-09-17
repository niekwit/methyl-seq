# Stage optional user-defined custom regions (config boxplot:regions) into
# the standard bed/{region}.bed location used by every region-wildcarded
# rule below. wildcard_constraints restricts this rule to exactly those
# custom names, so it can never collide with the standard region outputs
# from the generate_regions rule (resources.smk), which are fixed files
# rather than wildcard matches.
_custom_regions = config["boxplot"].get("regions", None) or {}

if _custom_regions:

    rule stage_custom_region_bed:
        input:
            lambda wildcards: _custom_regions[wildcards.region],
        output:
            "bed/{region}.bed",
        wildcard_constraints:
            region="|".join(re.escape(r) for r in _custom_regions),
        log:
            "logs/resources/stage_custom_region_{region}.log",
        threads: 1
        shell:
            "cp {input} {output} 2> {log}"


# Filter CpG probes to keep only those covered by >10 reads (per condition)
# -----------------------------------------------------
rule filter_cpg_probes_for_reads:
    input:
        cpg="resources/cpg_probes_sorted.bed",
        meth="results/bed/CpG_merged_{condition}.bed",
        cs="resources/chrom_sizes.txt",
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
        "bedtools intersect -wa -wb -sorted "
        "-a {input.cpg} -b {input.meth} -g {input.cs} | "
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
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        region_probes="resources/cpg_probes_{region}.bed",
    log:
        "logs/resources/filter_probes_{region}.log",
    threads: 1
    resources:
        runtime=15,
        mem_mb=2000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bedtools intersect -sorted -wa "
        "-a {input.probes} "
        "-b {input.regions} "
        "-g {input.chrom_sizes} | "
        "sort -k1,1 -k2,2n > {output.region_probes} 2> {log}"


# Create input for boxplot
# -----------------------------------------------------
rule boxplot_data:
    input:
        probes="resources/cpg_probes_{region}.bed",
        meth="results/bed/CpG_merged_{condition}.bed",
        cs="resources/chrom_sizes.txt",
    output:
        "results/boxplot/CpG_methylation_{condition}_{region}.txt",
    log:
        "logs/score_methylation_calls/boxplot_{condition}_{region}.log",
    threads: 1
    resources:
        runtime=360,
        mem_mb=2000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bedtools intersect -sorted -wa -wb -a {input.probes} -b {input.meth} -g {input.cs} | "
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
        csv="results/plots/boxplots_data.csv",
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


# Combine and plot LINE1 boxplot data (its own figure, separate from the
# main boxplots.pdf above) -- only defined when config boxplot:LINE1:plot
# is true. filter_cpg_probes_for_regions/boxplot_data above need no
# changes: they're already generic over {region}, and bed/LINE1.bed /
# bed/{subfamily}.bed (from generate_line1_regions in resources.smk) flow
# through them exactly like any other region, since LINE1_REGIONS is
# included in the Snakefile's global `region` wildcard_constraint.
# -----------------------------------------------------
if LINE1_REGIONS:

    rule combine_line1_boxplot_data:
        input:
            data=expand(
                "results/boxplot/CpG_methylation_{condition}_{region}.txt",
                condition=CONDITIONS,
                region=LINE1_REGIONS,
            ),
        output:
            "results/boxplot/CpG_methylation_all_conditions_line1_regions.txt",
        log:
            "logs/boxplot/combine_line1_data.log",
        threads: 1
        resources:
            runtime=10,
            mem_mb=2000,
        conda:
            "../envs/deeptools.yaml"
        shell:
            "cat {input.data} | "
            r"sed 's/^\s*//;s/\s/\t/g' > {output}"

    use rule plot_boxplot as plot_line1_boxplot with:
        input:
            "results/boxplot/CpG_methylation_all_conditions_line1_regions.txt",
        output:
            pdf="results/plots/line1_boxplots.pdf",
            csv="results/plots/line1_boxplots_data.csv",
        params:
            regions=LINE1_REGIONS,
        log:
            "logs/boxplot/plot_line1_boxplots.log",
