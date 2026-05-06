# Convert BedGraph to BigWig
# -----------------------------------------------------
rule bedgraph_to_bigwig:
    input:
        bg="results/bismark/{sample}/{sample}.deduplicated.bedGraph.gz",
        cs="resources/chrom_sizes.txt",
    output:
        bw="results/bismark/{sample}/{sample}.deduplicated.bw",
        bg=temp("results/temp/{sample}.bg"),
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


# Decompress bedGraph files for averaging
# -----------------------------------------------------
#rule decompress_bedgraph:
#    input:
#        bg="results/bismark/{sample}/{sample}.deduplicated.bedGraph.gz",
#    output:
#        # Change extension to .bg for wiggletools compatibility
#        bg=temp("results/temp/{sample}.bg"),
#    resources:
#        runtime=15,
#        mem_mb=4000,
#    threads: 4
#    log:
#        "logs/decompress_bedgraph/{sample}.log",
#    conda:
#        "../envs/deeptools.yaml"
#    shell:
#        "pigz -p {threads} -dc {input.bg} > {output.bg} 2> {log}"


# Create average BedGraph files
# -----------------------------------------------------
rule average_bedgraphs:
    input:
        # All replicate bedGraphs for the current condition
        bg=lambda wildcards: expand(
            "results/temp/{sample}.bg",
            sample=[s for s in SAMPLES if s.startswith(wildcards.condition)]
        ),
    output:
        bg=temp("results/temp/{condition}.bg"),
    resources:
        runtime=120,
    log:
        "logs/average_bedgraphs/{condition}.log",
    threads: 2
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/average_bedgraph.py"


# Sort average BedGraph files
# -----------------------------------------------------
rule sort_bedgraph:
    input:
        bg="results/temp/{condition}.bg",
    output:
        bg="results/bismark/{condition}.bg",
    log:
        "logs/sort_bedgraph/{condition}.log",
    threads: 1
    resources:
        runtime=30,
        mem_mb=2000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "LC_ALL=C sort -k1,1 -k2,2n {input.bg} > {output.bg}"


# Convert average BedGraph to BigWig
# -----------------------------------------------------
rule average_bedgraph_to_bigwig:
    input:
        bg="results/bismark/{condition}.bg",
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


# Create per-sample coverage bedGraph from Bismark coverage file
# -----------------------------------------------------
rule coverage_bedgraph:
    input:
        cov="results/bismark/{sample}/{sample}.deduplicated.bismark.cov.gz",
    output:
        bg=temp("results/bismark/{sample}/coverage/{sample}.bg"),
    log:
        "logs/coverage_bedgraph/{sample}.log",
    threads: 2
    resources:
        runtime=30,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "zcat {input.cov} | "
        """awk 'OFS="\\t" {{ print $1, $2, $2+1, $5+$6 }}' | """
        "LC_COLLATE=C sort -k1,1 -k2,2n > {output.bg} 2> {log}"


# Average coverage bedGraphs across replicates per condition
# -----------------------------------------------------
rule average_coverage_bedgraphs:
    input:
        bg=lambda wildcards: expand(
            "results/bismark/{sample}/coverage/{sample}.bg",
            sample=[s for s in SAMPLES if s.startswith(wildcards.condition)],
        ),
    output:
        bg=temp("results/bismark/coverage/{condition}.bg"),
    log:
        "logs/average_coverage_bedgraphs/{condition}.log",
    threads: 2
    resources:
        runtime=120,
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/average_bedgraph.py"


# Sort average coverage bedGraph for bedGraphToBigWig
# -----------------------------------------------------
rule sort_coverage_bedgraph:
    input:
        bg="results/bismark/coverage/{condition}.bg",
    output:
        bg=temp("results/bismark/coverage/sorted/{condition}.bg"),
    log:
        "logs/sort_coverage_bedgraph/{condition}.log",
    threads: 1
    resources:
        runtime=30,
        mem_mb=2000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "LC_ALL=C sort -k1,1 -k2,2n {input.bg} > {output.bg} 2> {log}"


# Convert average coverage bedGraph to BigWig
# -----------------------------------------------------
rule average_coverage_bedgraph_to_bigwig:
    input:
        bg="results/bismark/coverage/sorted/{condition}.bg",
        cs="resources/chrom_sizes.txt",
    output:
        "results/bigwig/coverage/{condition}.bw",
    log:
        "logs/average_coverage_bedgraph_to_bigwig/{condition}.log",
    threads: 2
    resources:
        runtime=60,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bedGraphToBigWig {input.bg} {input.cs} {output} 2> {log}"