# CpG methylation profiles/heatmaps (deepTools computeMatrix / plotHeatmap
# / plotProfile) over results/bigwig/{condition}.bw, for two kinds of
# regions:
#   - DMRs (reference-point, centred on each DMR's midpoint), split into
#     separate hypo-/hypermethylated runs.
#   - TE boxplot regions (scale-regions, --regionBodyLength = the class
#     block's own min_length), one combined plot per configured TE class
#     block (class total + configured family/subfamilies together).
#
# Neither set of rules is wrapped in an `if config["DMR"]["run"]`/
# `if TE_REGIONS` guard: targets() (common.smk) already only requests
# these outputs when the corresponding feature is enabled, and -- unlike
# e.g. resources.smk's `if _te_active` guard, which exists to stop
# subtract_te_from_genic/... from colliding with generate_regions' own
# output path for the same filename when TE is inactive -- none of these
# rules' output paths are ever produced by another rule, so there's no
# ambiguity to avoid by hiding the rule definitions themselves. The TE
# per-class rules below are also naturally a no-op (the for loop just
# doesn't iterate) when no TE class block is configured, without needing
# a guard of their own.


# results/dmrs/{status}_DMRs.bed -- see identify_dmrs (dmrs.smk)
rule dmr_profile_matrix:
    input:
        bw=expand("results/bigwig/{condition}.bw", condition=CONDITIONS),
        bed="results/dmrs/{status}_DMRs.bed",
    output:
        "results/deeptools/dmr_{status}_matrix.gz",
    wildcard_constraints:
        status="hypomethylated|hypermethylated",
    params:
        upstream=config["DMR"]["profile"]["upstream"],
        downstream=config["DMR"]["profile"]["downstream"],
        binSize=config["DMR"]["profile"]["binSize"],
        extra=config["DMR"]["profile"]["extra"],
    log:
        "logs/deeptools/dmr_profile_matrix_{status}.log",
    threads: 8
    resources:
        runtime=60,
        mem_mb=8000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "computeMatrix reference-point "
        "--referencePoint center "
        "-S {input.bw} "
        "-R {input.bed} "
        "-a {params.upstream} "
        "-b {params.downstream} "
        "--skipZeros "
        "--binSize {params.binSize} "
        "--averageTypeBins mean "
        "--numberOfProcessors {threads} "
        "{params.extra} "
        "--outFileName {output} "
        "2> {log}"


rule dmr_profile_heatmap:
    input:
        "results/deeptools/dmr_{status}_matrix.gz",
    output:
        "results/plots/dmrs/{status}_heatmap.pdf",
    wildcard_constraints:
        status="hypomethylated|hypermethylated",
    params:
        samples_label=CONDITIONS,
    log:
        "logs/deeptools/dmr_heatmap_{status}.log",
    threads: 1
    resources:
        runtime=15,
        mem_mb=4000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotHeatmap "
        "--matrixFile {input} "
        "--outFileName {output} "
        "--heatmapHeight 10 "
        "--heatmapWidth 5 "
        "--perGroup "
        "--samplesLabel {params.samples_label} "
        "--regionsLabel '{wildcards.status} DMRs' "
        "--yAxisLabel 'CpG methylation (%)' "
        "--averageType mean "
        "--plotTitle '{wildcards.status} DMRs' "
        "--colorMap viridis "
        "2> {log}"


rule dmr_profile_plot:
    input:
        "results/deeptools/dmr_{status}_matrix.gz",
    output:
        "results/plots/dmrs/{status}_profile.pdf",
    wildcard_constraints:
        status="hypomethylated|hypermethylated",
    params:
        samples_label=CONDITIONS,
    log:
        "logs/deeptools/dmr_profile_{status}.log",
    threads: 1
    resources:
        runtime=15,
        mem_mb=4000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotProfile "
        "--matrixFile {input} "
        "--outFileName {output} "
        "--perGroup "
        "--samplesLabel {params.samples_label} "
        "--regionsLabel '{wildcards.status} DMRs' "
        "--yAxisLabel 'CpG methylation (%)' "
        "--averageType mean "
        "--plotTitle '{wildcards.status} DMRs' "
        "2> {log}"


for _class_name, _block in te_class_blocks():
    _region_names = te_class_region_list(_class_name, _block)
    _profile = te_profile_settings(_block)

    rule:
        name:
            f"te_profile_matrix_{_class_name}"
        input:
            bw=expand("results/bigwig/{condition}.bw", condition=CONDITIONS),
            bed=expand("bed/{region}.bed", region=_region_names),
        output:
            f"results/deeptools/te_{_class_name}_matrix.gz",
        params:
            region_body_length=_block["min_length"],
            upstream=_profile["upstream"],
            downstream=_profile["downstream"],
            binSize=_profile["binSize"],
            extra=_profile["extra"],
        log:
            f"logs/deeptools/te_profile_matrix_{_class_name}.log",
        threads: 8
        resources:
            runtime=60,
            mem_mb=8000,
        conda:
            "../envs/deeptools.yaml"
        shell:
            "computeMatrix scale-regions "
            "-S {input.bw} "
            "-R {input.bed} "
            "-a {params.upstream} "
            "-b {params.downstream} "
            "--regionBodyLength {params.region_body_length} "
            "--smartLabels "
            "--skipZeros "
            "--binSize {params.binSize} "
            "--averageTypeBins mean "
            "--numberOfProcessors {threads} "
            "{params.extra} "
            "--outFileName {output} "
            "2> {log}"

    rule:
        name:
            f"te_profile_heatmap_{_class_name}"
        input:
            f"results/deeptools/te_{_class_name}_matrix.gz",
        output:
            f"results/plots/te_{_class_name}_heatmap.pdf",
        params:
            samples_label=CONDITIONS,
            regions_label=_region_names,
        log:
            f"logs/deeptools/te_heatmap_{_class_name}.log",
        threads: 1
        resources:
            runtime=15,
            mem_mb=4000,
        conda:
            "../envs/deeptools.yaml"
        shell:
            "plotHeatmap "
            "--matrixFile {input} "
            "--outFileName {output} "
            "--heatmapHeight 10 "
            "--heatmapWidth 5 "
            "--perGroup "
            "--samplesLabel {params.samples_label} "
            "--regionsLabel {params.regions_label} "
            "--yAxisLabel 'CpG methylation (%)' "
            "--averageType mean "
            f"--plotTitle 'CpG methylation at {_class_name} elements' "
            "--colorMap viridis "
            "2> {log}"
