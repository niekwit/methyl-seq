# Subset average bedGraph for all CpGs that overlap with
# specified genomic regions (i.e. BED files)
# -----------------------------------------------------
rule metaplot_data:
    input:
        bg="results/bismark/{condition}.bg",
        regions="bed/{meta_region}_extended.bed",
    output:
        "results/metaplot/{condition}_{meta_region}.bg",
    log:
        "logs/metaplot/data_{condition}_{meta_region}.log",
    threads: 1
    resources:
        runtime=180,
        mem_mb=2000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bedtools intersect -wa -wb "
        "-a {input.bg} "
        "-b {input.regions} | "
        "sort -k1,1 -k2,2n "
        "sed 's/$/\t{wildcards.condition}/' > {output} 2> {log}"


# Create metaplot using R
# -----------------------------------------------------
rule plot_metaplot:
    input:
        bg=expand(
            "results/metaplot/{condition}_{meta_region}.bg",
            condition=CONDITIONS,
            meta_region=META_REGIONS,
        ),
    output:
        "results/plots/metaplots.pdf",
    params:
        extend_with=config["metaplot"]["extend"],
    log:
        "logs/metaplot/plot.log",
    threads: 4
    resources:
        runtime=30,
        mem_mb=8000,
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_metaplot.R"
