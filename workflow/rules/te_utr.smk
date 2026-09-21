# Optional LINE1 5' UTR vs. remainder %CpG methylation boxplot, opt-in per
# TE class block via a `utr_analysis:` sub-block (only valid when that
# block's family is exactly "L1" -- see validate_te_utr_config() in
# common.smk). Each near-full-length L1 element is split into a 5' UTR and
# a "remainder" (ORF1 + linker + ORF2 + 3' UTR) by locating ORF1/ORF2 via
# minimap2 against a user-supplied GenBank consensus record
# (line1_annotate_utrs_orfs.py), and each element's UTR/remainder span
# becomes its own boxplot data point -- bypassing the CpG-probe pooling
# layer used by the standard region/TE boxplots, since a probe-sized chunk
# would be a biologically meaningless split of one element's UTR.
_UTR_CLASS_BLOCKS = te_utr_class_blocks()

if _UTR_CLASS_BLOCKS:

    for _class_name, _block in _UTR_CLASS_BLOCKS:
        _utr = _block["utr_analysis"]
        _family_name = te_family_list(_block)[0]  # == "L1", validated
        _te_names = [_family_name] + list(_block.get("subfamilies", []))

        rule:
            name:
                f"annotate_line1_utrs_{_class_name}"
            input:
                bed=f"bed/{_family_name}.bed",
                genbank=_utr["genbank"],
                fasta="resources/combined_genome.fa",
                fai="resources/combined_genome.fa.fai",
            output:
                f"bed/{_class_name}_UTRs_ORFs.bed",
            params:
                script=os.path.join(WORKFLOW_SCRIPTS, "line1_annotate_utrs_orfs.py"),
                min_mapq=_utr.get("min_mapq", 20),
                min_identity=_utr.get("min_identity", 0.6),
                min_coverage=_utr.get("min_coverage", 0.5),
                single_orf_flag=" --single-orf" if _utr.get("single_orf", False) else "",
            log:
                f"logs/resources/annotate_line1_utrs_{_class_name}.log",
            conda:
                "../envs/line1_utr.yaml"
            threads: 1
            resources:
                runtime=360,
                mem_mb=8000,
            shell:
                "python {params.script} "
                "--bed {input.bed} "
                "--genbank {input.genbank} "
                "--fasta {input.fasta} "
                "--out {output} "
                "--min-length 0 "
                "--min-mapq {params.min_mapq} "
                "--min-identity {params.min_identity} "
                "--min-coverage {params.min_coverage} "
                "{params.single_orf_flag} "
                "--log {log}"

        for _te_name in _te_names:
            _source_bed = f"bed/{_te_name}.bed"

            rule:
                name:
                    f"generate_line1_utr_regions_{_te_name}"
                input:
                    annotated=f"bed/{_class_name}_UTRs_ORFs.bed",
                    te=_source_bed,
                output:
                    f"bed/{_te_name}_5UTR_remainder.bed",
                log:
                    f"logs/resources/generate_line1_utr_regions_{_te_name}.log",
                conda:
                    "../envs/deeptools.yaml"
                threads: 1
                resources:
                    runtime=15,
                    mem_mb=2000,
                shell:
                    "bedtools intersect -u -a {input.annotated} -b {input.te} | "
                    r"awk '$4 ~ /_(5UTR|remainder)$/' | "
                    "sort -k1,1 -k2,2n "
                    "> {output} 2> {log}"

    # Per (te_name, condition) methylation-call counting -- mirrors the
    # reference box_plot_data_te_5utr_remainder.sh exactly: intersect each
    # element's UTR/remainder span directly against the per-condition
    # merged CpG calls (results/bed/CpG_merged_{condition}.bed, 6 columns:
    # chrom, start, end, read ID, methylation call Z/z, strand), so column
    # 4 of the -wa/-wb output is the UTR/remainder feature name and column
    # 11 (6 + the meth file's own column 5) is the methylation call.
    rule line1_utr_boxplot_data:
        input:
            bed="bed/{te_name}_5UTR_remainder.bed",
            meth="results/bed/CpG_merged_{condition}.bed",
            cs="resources/chrom_sizes.txt",
        output:
            "results/boxplot/te_5utr/{te_name}_{condition}.txt",
        wildcard_constraints:
            te_name="|".join(re.escape(t) for t in UTR_TE_NAMES),
        log:
            "logs/boxplot/te_5utr_{te_name}_{condition}.log",
        threads: 1
        resources:
            runtime=60,
            mem_mb=2000,
        conda:
            "../envs/deeptools.yaml"
        shell:
            "bedtools intersect -sorted -wa -wb -a {input.bed} -b {input.meth} -g {input.cs} | "
            "cut -f4,11 | "
            "sort -k1,1 -k2,2n | "
            "uniq -c | "
            "sed 's/$/\t{wildcards.te_name}\t{wildcards.condition}/' > {output}"

    rule combine_line1_utr_boxplot_data:
        input:
            expand(
                "results/boxplot/te_5utr/{te_name}_{condition}.txt",
                te_name=UTR_TE_NAMES,
                condition=CONDITIONS,
            ),
        output:
            "results/boxplot/CpG_methylation_all_conditions_te_5utr.txt",
        log:
            "logs/boxplot/combine_te_5utr_data.log",
        threads: 1
        resources:
            runtime=10,
            mem_mb=2000,
        conda:
            "../envs/deeptools.yaml"
        shell:
            "cat {input} | "
            r"sed 's/^\s*//;s/\s/\t/g' > {output}"

    rule plot_line1_utr_boxplot:
        input:
            "results/boxplot/CpG_methylation_all_conditions_te_5utr.txt",
        output:
            pdf="results/plots/te_5utr_boxplots.pdf",
            csv="results/plots/te_5utr_boxplots_data.csv",
        params:
            te_names=UTR_TE_NAMES,
            min_reads=config["boxplot"]["min_reads"],
        log:
            "logs/boxplot/plot_te_5utr_boxplots.log",
        threads: 1
        resources:
            runtime=30,
            mem_mb=2000,
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/plot_te_5utr_boxplot.R"

    # 5' UTR methylation-status classification: mean %CpG methylation over
    # just each element's 5' UTR span, per condition, taken from the
    # per-condition averaged bigwig (results/bigwig/{condition}.bw -- see
    # average_bedgraph_to_bigwig in bigwig.smk) -- mirrors icr_scores in
    # icr.smk, restricted to the "_5UTR" rows of each TE's own
    # bed/{te_name}_5UTR_remainder.bed ("_remainder" rows are irrelevant
    # to this classification).
    rule line1_5utr_bigwig_average:
        input:
            bed="bed/{te_name}_5UTR_remainder.bed",
            bw="results/bigwig/{condition}.bw",
        output:
            temp("results/te_5utr_status/{te_name}_{condition}_avg.txt"),
        wildcard_constraints:
            te_name="|".join(re.escape(t) for t in UTR_TE_NAMES),
        log:
            "logs/te_5utr_status/{te_name}_{condition}_avg.log",
        threads: 1
        resources:
            runtime=15,
            mem_mb=2000,
        conda:
            "../envs/deeptools.yaml"
        shell:
            # bigWigAverageOverBed doesn't accept the trailing score/strand
            # columns (score="." isn't a valid integer to it) -- cut down to
            # plain chrom/start/end/name first, as the reference script did.
            r"awk '$4 ~ /_5UTR$/' {input.bed} | "
            "cut -f1-4 | "
            "bigWigAverageOverBed {input.bw} /dev/stdin /dev/stdout 2> {log} | "
            r"sed 's/$/\t{wildcards.te_name}\t{wildcards.condition}/' > {output}"

    rule combine_line1_5utr_bigwig_average:
        input:
            expand(
                "results/te_5utr_status/{te_name}_{condition}_avg.txt",
                te_name=UTR_TE_NAMES,
                condition=CONDITIONS,
            ),
        output:
            "results/te_5utr_status/all_conditions_5utr_avg.txt",
        log:
            "logs/te_5utr_status/combine_5utr_avg.log",
        threads: 1
        resources:
            runtime=10,
            mem_mb=2000,
        conda:
            "../envs/deeptools.yaml"
        shell:
            "cat {input} > {output} 2> {log}"

    rule plot_line1_5utr_methylation_status:
        input:
            avg="results/te_5utr_status/all_conditions_5utr_avg.txt",
            beds=expand("bed/{te_name}_5UTR_remainder.bed", te_name=UTR_TE_NAMES),
        output:
            pdf="results/plots/te_5utr_methylation_histogram.pdf",
            csv="results/plots/te_5utr_methylation_histogram_data.csv",
            # Cross-condition classification, one comparison per non-reference
            # ("KO") condition against the single config DMR:reference_condition
            # ("REF"): active = REF >= cutoff & KO < cutoff (lost methylation
            # upon KO), inactive = REF >= cutoff & KO >= cutoff (stays
            # methylated); elements with REF < cutoff are excluded from both.
            active=expand(
                "results/te_5utr_status/{te_name}_{ko_condition}_active.bed",
                te_name=UTR_TE_NAMES,
                ko_condition=UTR_TE_KO_CONDITIONS,
            ),
            inactive=expand(
                "results/te_5utr_status/{te_name}_{ko_condition}_inactive.bed",
                te_name=UTR_TE_NAMES,
                ko_condition=UTR_TE_KO_CONDITIONS,
            ),
        params:
            te_names=UTR_TE_NAMES,
            reference_condition=config["DMR"]["reference_condition"],
            ko_conditions=UTR_TE_KO_CONDITIONS,
            cutoffs=[UTR_TE_CUTOFFS[t] for t in UTR_TE_NAMES],
        log:
            "logs/te_5utr_status/plot_methylation_status.log",
        threads: 1
        resources:
            runtime=30,
            mem_mb=2000,
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/plot_te_5utr_methylation_status.R"
