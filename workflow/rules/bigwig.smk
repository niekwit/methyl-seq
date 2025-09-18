# Convert BedGraph to BigWig
# -----------------------------------------------------
rule bedgraph_to_bigwig:
    input:
        bg="results/bismark/{sample}/{sample}.deduplicated.bedGraph.gz",
        cs="resources/chrom_sizes.txt",
    output:
        bw="results/bismark/{sample}/{sample}.deduplicated.bw",
        bg=temp("results/bismark/{sample}/{sample}.deduplicated.bedGraph"),
    log:
        "logs/bedgraph_to_bigwig/{sample}.log",
    threads: 4
    resources:
        runtime=120,
    conda:
        "../envs/deeptools.yaml"
    shell:
        # bedGraphToBigWig does not like piped input
        "pigz -p 4 -dc {input.bg} > {output.bg}; "
        "bedGraphToBigWig "
        "{output.bg} "
        "{input.cs} "
        "{output.bw} "
        "2> {log}"


# Create summary of all BigWig files
# -----------------------------------------------------
rule bigwig_summary:
    input:
        bw=expand("results/bismark/{sample}/{sample}.deduplicated.bw", sample=SAMPLES),
    output:
        "results/deeptools/bigwig_summary.npz",
    log:
        "logs/deeptools/bigwig_summary.log",
    threads: 12
    resources:
        runtime=120,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "multiBigwigSummary bins "
        "--bwfiles {input.bw} "
        "--outFile {output} "
        "2> {log}"


# PCA on BigWig files
# -----------------------------------------------------
rule PCA:
    input:
        "results/deeptools/bigwig_summary.npz",
    output:
        "results/deeptools/PCA.tab",
    params:
        extra=config["deeptools"]["plotPCA"]["extra"],
    threads: 8
    resources:
        runtime=60,
    log:
        "logs/deeptools/PCA.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotPCA "
        "--corData {input} "
        "--outFileNameData {output} "
        "--transpose "
        "{params.extra} "
        "> {log} 2>&1"


# Plot PCA results
# -----------------------------------------------------
rule plotPCA:
    input:
        "results/deeptools/PCA.tab",
    output:
        pca=report("results/plots/PCA.pdf", caption="../report/pca.rst", category="PCA"),
        scree=report(
            "results/plots/scree.pdf", caption="../report/scree.rst", category="PCA"
        ),
    params:
        extra="",
    threads: 1
    resources:
        runtime=15,
    log:
        "logs/plotting/plotPCA.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_PCA.R"


# Create average BedGraph files
# -----------------------------------------------------
rule average_bedgraphs:
    input:
        # All replicate bedGraphs for the current condition
        bg=lambda wildcards: expand(
            "results/bismark/{sample}/{sample}.deduplicated.bedGraph.gz",
            sample=[s for s in SAMPLES if s.startswith(wildcards.condition)]
        )
    output:
        bg=temp("results/bismark/{condition}_average_deduplicated.bedGraph"),
    resources:
        runtime=120,
    log:
        "logs/average_bedgraphs/{condition}.log",
    threads: 2
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/average_bedgraph.py"


# Convert average BedGraph to BigWig
# -----------------------------------------------------
rule average_bedgraph_to_bigwig:
    input:
        bg="results/bismark/{condition}_average_deduplicated.bedGraph",
        cs="resources/chrom_sizes.txt",
    output:
        "results/bigwig/{condition}.bw",
    threads: 2
    resources:
        runtime=60,
    log:
        "logs/average_bedgraph_to_bigwig/{condition}.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bedGraphToBigWig {input.bg} {input.cs} {output} 2> {log}"
