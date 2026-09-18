# Paternally imprinted region (ICR) methylation heatmap. mm39 only for now
# (resources.icr_regions is None for every other genome -- see
# resources.py) -- the "test" genome gets a single-region (Rasgrf1 only)
# subset so this code path is exercised in CI.
if resources.icr_regions:

    if resources.icr_regions_url:

        # Only the "test" genome needs this: its ICR bed is a small,
        # already-checked-in fixture fetched the same way as every other
        # test-genome resource (see resources.py). mm39's icr_regions
        # points directly at the static, checked-in
        # workflow/resources/icr_regions_mm39.bed instead -- no download
        # rule needed for that one.
        rule get_icr_regions:
            output:
                resources.icr_regions,
            params:
                url=resources.icr_regions_url,
            log:
                "logs/resources/get_icr_regions.log",
            conda:
                "../envs/bismark.yaml"
            threads: 1
            shell:
                "wget -q {params.url} -O {output} 2> {log}"

    # Average %CpG methylation per ICR, per condition
    rule icr_scores:
        input:
            bw="results/bigwig/{condition}.bw",
            icr=resources.icr_regions,
        output:
            temp("results/icr/{condition}_scores.txt"),
        log:
            "logs/icr/{condition}_scores.log",
        threads: 1
        resources:
            runtime=15,
            mem_mb=2000,
        conda:
            "../envs/deeptools.yaml"
        shell:
            "bigWigAverageOverBed {input.bw} {input.icr} /dev/stdout 2> {log} | "
            "sed 's/$/\t{wildcards.condition}/' > {output}"

    rule combine_icr_scores:
        input:
            expand("results/icr/{condition}_scores.txt", condition=CONDITIONS),
        output:
            "results/icr/all_conditions_scores.txt",
        log:
            "logs/icr/combine_scores.log",
        threads: 1
        resources:
            runtime=10,
            mem_mb=2000,
        conda:
            "../envs/deeptools.yaml"
        shell:
            "cat {input} > {output} 2> {log}"

    rule plot_icr_heatmap:
        input:
            scores="results/icr/all_conditions_scores.txt",
            icr_bed=resources.icr_regions,
        output:
            pdf="results/plots/icr_heatmap.pdf",
            csv="results/plots/icr_heatmap_data.csv",
        log:
            "logs/icr/plot_icr_heatmap.log",
        threads: 1
        resources:
            runtime=15,
            mem_mb=2000,
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/icr_heatmap.R"
