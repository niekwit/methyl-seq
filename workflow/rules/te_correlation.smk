# Sample-level CpG methylation correlation over the elements of each
# configured TE family (every `family` value under any TE class block in
# config boxplot, e.g. L1): deepTools multiBigwigSummary BED-file over the
# family's own bed/{family}.bed (the same min_length- and gene-body-filtered
# element set as the family's boxplot/heatmap), followed by a Spearman
# plotCorrelation heatmap. Uses the per-sample bigWigs (not the
# per-condition averages) so replicates can be compared with each other.
#
# Nothing is defined when no TE class block configures a `family`
# (TE_FAMILIES == []): targets() (common.smk) only requests these outputs
# for the configured families, and the family wildcard_constraint below
# can't be built from an empty list.
if TE_FAMILIES:

    rule te_family_bigwig_summary:
        input:
            bw=expand(
                "results/bismark/{sample}/{sample}.deduplicated.bw", sample=SAMPLES
            ),
            bed="bed/{family}.bed",
        output:
            "results/deeptools/te_{family}_bigwig_summary.npz",
        wildcard_constraints:
            family="|".join(re.escape(f) for f in TE_FAMILIES),
        params:
            labels=SAMPLES,
        log:
            "logs/deeptools/te_{family}_bigwig_summary.log",
        threads: 4
        resources:
            runtime=60,
            mem_mb=8000,
        conda:
            "../envs/deeptools.yaml"
        shell:
            "multiBigwigSummary BED-file "
            "--bwfiles {input.bw} "
            "--BED {input.bed} "
            "--labels {params.labels} "
            "--numberOfProcessors {threads} "
            "--outFileName {output} "
            "2> {log}"

    rule te_family_correlation:
        input:
            "results/deeptools/te_{family}_bigwig_summary.npz",
        output:
            pdf="results/plots/te_{family}_correlation.pdf",
            matrix="results/deeptools/te_{family}_correlation.tab",
        wildcard_constraints:
            family="|".join(re.escape(f) for f in TE_FAMILIES),
        log:
            "logs/deeptools/te_{family}_correlation.log",
        threads: 1
        resources:
            runtime=15,
            mem_mb=4000,
        conda:
            "../envs/deeptools.yaml"
        shell:
            "plotCorrelation "
            "--corData {input} "
            "--corMethod spearman "
            "--whatToPlot heatmap "
            "--removeOutliers "
            "--plotNumbers "
            "--colorMap viridis "
            "--plotTitle 'Spearman correlation of CpG methylation at {wildcards.family} elements' "
            "--outFileCorMatrix {output.matrix} "
            "--plotFile {output.pdf} "
            "2> {log}"
