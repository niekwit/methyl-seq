# Identify and annotate DMRs using methylKit
# -----------------------------------------------------
rule identify_dmrs:
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
        hypo="results/dmrs/hypomethylated_DMRs.bed",
        hyper="results/dmrs/hypermethylated_DMRs.bed",
        meth_tiles="results/dmrs/all_methylation_tiles.rds",
        diff_tiles="results/dmrs/differential_methylation_tiles.rds",
        sign_diff_tiles="results/dmrs/significant_differential_methylation_tiles.rds",
    params:
        samples=SAMPLES,
        tile_size=config["DMR"]["tile_size"],
        step_size=config["DMR"]["step_size"],
        min_per_group=config["DMR"]["min_per_group"],
        ref_cond=config["DMR"]["reference_condition"],
    log:
        "logs/dmrs/identify_dmrs.log",
    threads: 14
    resources:
        runtime=120,
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/dmr.R"


# Annotate DMRs using ChIPseeker
# -----------------------------------------------------
rule annotate_dmrs:
    input:
        hypo="results/dmrs/hypomethylated_DMRs.bed",
        hyper="results/dmrs/hypermethylated_DMRs.bed",
        diff_tiles="results/dmrs/differential_methylation_tiles.rds",
    output:
        hypo="results/dmrs/hypomethylated_DMRs_annotated.tab",
        hyper="results/dmrs/hypermethylated_DMRs_annotated.tab",
        distance="results/plots/dmrs/DMR_distance_to_TSS.pdf",
        distribution="results/plots/dmrs/DMR_genomic_distribution.pdf",
        volcano="results/plots/dmrs/DMR_volcano.pdf",
    log:
        "logs/dmrs/annotate_dmrs.log",
    threads: 4
    resources:
        runtime=60,
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/annotate_dmrs.R"
