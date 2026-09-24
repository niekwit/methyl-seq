# Absolute paths into the workflow's own scripts/resources -- shell commands
# run with CWD = --directory (an analysis dir, not necessarily this repo),
# so a bare relative "workflow/scripts/foo.py" only resolves by accident.
WORKFLOW_SCRIPTS = os.path.join(workflow.basedir, "scripts")
WORKFLOW_RESOURCES = os.path.join(workflow.basedir, "resources")


rule get_genome_fasta:
    output:
        resources.fasta,
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log",
    conda:
        "../envs/bismark.yaml"
    threads: 1
    shell:
        "wget -q {params.url} -O {output}.gz 2> {log} ;"
        "pigz -d -c {output}.gz > {output}"


rule filter_fasta:
    input:
        fastx=resources.fasta,
    output:
        fastx=resources.filtered_fasta,
    log:
        "logs/seqkit/filter_fasta.log",
    params:
        command="grep",
        extra="-r -p '^([0-9]+|X|Y|MT)$'",  # keep only standard chromosomes
    threads: 4
    wrapper:
        "v9.4.2/bio/seqkit"


rule get_control_fasta:
    output:
        resources.control_fasta,
    params:
        url=resources.control_fasta_url,
    retries: 3
    log:
        "logs/resources/get_control_fasta.log",
    conda:
        "../envs/bismark.yaml"
    threads: 1
    shell:
        "wget -q {params.url} -O {output} 2> {log}"


rule combine_fasta:
    input:
        genome=resources.filtered_fasta,
        control=resources.control_fasta,
    output:
        "resources/combined_genome.fa",
    log:
        "logs/resources/combine_fasta.log",
    threads: 1
    conda:
        "../envs/bismark.yaml"
    shell:
        "cat {input.genome} {input.control}  > {output} 2> {log}"


rule index_fasta:
    input:
        "resources/combined_genome.fa",
    output:
        "resources/combined_genome.fa.fai",
    log:
        "logs/resources/index_fasta.log",
    threads: 1
    conda:
        "../envs/bismark.yaml"
    shell:
        "samtools faidx {input} 2> {log}"


rule chrom_sizes:
    input:
        "resources/combined_genome.fa.fai",
    output:
        "resources/chrom_sizes.txt",
    log:
        "logs/resources/chrom_sizes.log",
    threads: 1
    conda:
        "../envs/bismark.yaml"
    shell:
        "cut -f1,2 {input} > {output} 2> {log}"


# EM-seq spike-in/conversion-control contigs (methylation_controls.fa,
# concatenated onto the real genome in combine_fasta above), as a BED --
# same contig names already excluded from bed/whole_genome.bed
# (generate_regions.py), results/bed/CpG_merged_{condition}.bed
# (process_methylation_calls.smk), and dmr.R. Used to blacklist them from
# bigwig_summary's genome-wide binning (bigwig.smk), which otherwise has
# no per-region restriction to exclude them via.
rule control_dna_regions:
    input:
        "resources/chrom_sizes.txt",
    output:
        "resources/control_dna_regions.bed",
    log:
        "logs/resources/control_dna_regions.log",
    threads: 1
    shell:
        "awk -v OFS='\t' "
        "'$1 ~ /^(phage_T4|phage_Xp12|phage_lambda|plasmid_puc19c)$/ "
        "{{ print $1, 0, $2 }}' "
        "{input} > {output} 2> {log}"


use rule get_genome_fasta as get_gtf with:
    output:
        resources.gtf,
    params:
        url=resources.gtf_url,
    log:
        "logs/resources/get_gtf.log",


# regulatory_gtf/cpg_islands/repeat_mask are None for genomes with no such
# source (e.g. dm6, or the "test" genome's mm39 mini locus -- see
# resources.py), in which case these rules are not even defined: a rule
# can't declare `output: None`, and nothing downstream needs them (see the
# `if resources.regulatory_gtf:`-style guards in generate_regions below and
# in regions() in common.smk).
if resources.regulatory_gtf_url:

    use rule get_genome_fasta as get_regulatory_gtf with:
        output:
            resources.regulatory_gtf,
        params:
            url=resources.regulatory_gtf_url,
        log:
            "logs/resources/get_regulatory_gtf.log",


if resources.cpg_islands_url:

    # Downloads the UCSC cpgIslandExt table (goldenPath database dump -- the
    # same table hgTables serves under Track: CpG Islands) and converts it
    # to a plain BED file, stripping the "chr" prefix and keeping only
    # standard chromosomes to match the Ensembl-style contigs used
    # elsewhere in this workflow.
    rule prepare_cpg_islands:
        output:
            resources.cpg_islands,
        params:
            url=resources.cpg_islands_url,
        log:
            "logs/resources/prepare_cpg_islands.log",
        conda:
            "../envs/bismark.yaml"
        threads: 1
        shell:
            "wget -q {params.url} -O {output}.gz 2> {log} ;"
            "zcat {output}.gz | "
            "awk -F'\\t' -v OFS='\\t' "
            "'{{sub(/^chr/, \"\", $2); if ($2 ~ /^([0-9]+|X|Y|MT)$/) print $2, $3, $4, $5}}' "
            "> {output} 2>> {log}"


if resources.repeat_mask_url:

    # Downloads the UCSC rmsk table (goldenPath database dump -- the same
    # table hgTables serves under Track: RepeatMasker), converts it to a
    # plain BED file (chrom, start, end, repName, ., strand, repClass,
    # repFamily), stripping the "chr" prefix and keeping only standard
    # chromosomes to match the Ensembl-style contigs used elsewhere in this
    # workflow, and drops rows whose repClass (rmsk column 12) is a
    # non-transposable-element class listed in
    # workflow/resources/nonTE_repClasses.txt. Filtering is done on the
    # repClass field specifically (not a whole-line substring match), since
    # some genuinely transposable SINE families are themselves named after
    # non-TE classes (e.g. SINE/tRNA-RTE, SINE/tRNA-Deu are TEs, not tRNA
    # genes) and would otherwise be dropped by mistake.
    rule prepare_repeat_mask:
        input:
            nonte=os.path.join(WORKFLOW_RESOURCES, "nonTE_repClasses.txt"),
        output:
            resources.repeat_mask,
        params:
            url=resources.repeat_mask_url,
        log:
            "logs/resources/prepare_repeat_mask.log",
        conda:
            "../envs/bismark.yaml"
        threads: 1
        shell:
            "wget -q {params.url} -O {output}.gz 2> {log} ;"
            "zcat {output}.gz | "
            "awk -F'\\t' -v OFS='\\t' "
            "'NR==FNR {{nonte[$1]=1; next}} "
            '{{sub(/^chr/, "", $6); if ($6 !~ /^([0-9]+|X|Y|MT)$/) next; '
            "if ($12 in nonte) next; "
            'print $6, $7, $8, $11, ".", $10, $12, $13}}\' '
            "{input.nonte} - "
            "> {output} 2>> {log}"


# Standard genomic regions for boxplot/DMR annotation, generated via
# generate_regions.py. promoter/cpg_islands outputs are only declared when
# the genome actually has a regulatory GTF / CpG island track (see
# resources.py) -- see also regions() in common.smk, which builds REGIONS
# from these outputs plus any optional config-defined custom regions.
#
# When the TE boxplot feature is on (validate_te_config() in common.smk
# already guarantees resources.repeat_mask is set whenever TE_REGIONS is
# non-empty), genic/promoter/cpg_islands/intergenic are written to intermediate
# "_raw" paths here instead of their final bed/ location: the TE
# pipeline's gene-overlap exclusion (filter_repeat_mask_nongenic, below)
# needs the ORIGINAL, un-subtracted gene-body span, and the
# subtract_te_from_* rules below then produce the final, transposon-
# subtracted bed/{region}.bed that the existing boxplot pipeline consumes
# -- see the plan's "circularity to avoid" note for why these can't be the
# same file.
_te_active = bool(TE_REGIONS)

_region_outputs = {
    "whole_genome": "bed/whole_genome.bed",
    "genic": "resources/genic_raw.bed" if _te_active else "bed/genic.bed",
    "exon": "bed/exon.bed",
    "intron": "bed/intron.bed",
    "intergenic": (
        "resources/intergenic_raw.bed" if _te_active else "bed/intergenic.bed"
    ),
}
if resources.regulatory_gtf:
    _region_outputs["promoter"] = (
        "resources/promoter_raw.bed" if _te_active else "bed/promoter.bed"
    )
if resources.cpg_islands:
    _region_outputs["cpg_islands"] = (
        "resources/cpg_islands_raw.bed" if _te_active else "bed/cpg_islands.bed"
    )

_region_promoter_args = (
    f"--regulatory-gtf {resources.regulatory_gtf} --promoter-bed {_region_outputs['promoter']}"
    if resources.regulatory_gtf
    else ""
)
_region_cpg_args = (
    f"--cpg-island-bed {resources.cpg_islands} --cpg-island-bed-out {_region_outputs['cpg_islands']}"
    if resources.cpg_islands
    else ""
)


rule generate_regions:
    input:
        gtf=resources.gtf,
        chrom_sizes="resources/chrom_sizes.txt",
        regulatory_gtf=resources.regulatory_gtf if resources.regulatory_gtf else [],
        cpg_islands=resources.cpg_islands if resources.cpg_islands else [],
    output:
        **_region_outputs,
    params:
        script=os.path.join(WORKFLOW_SCRIPTS, "generate_regions.py"),
        promoter_args=_region_promoter_args,
        cpg_args=_region_cpg_args,
    log:
        "logs/resources/generate_regions.log",
    conda:
        "../envs/deeptools.yaml"
    threads: 1
    resources:
        runtime=30,
        mem_mb=8000,
    shell:
        "python {params.script} "
        "--gtf {input.gtf} "
        "--chrom-sizes {input.chrom_sizes} "
        "--whole-genome-bed {output.whole_genome} "
        "--genic-bed {output.genic} "
        "--exon-bed {output.exon} "
        "--intron-bed {output.intron} "
        "--intergenic-bed {output.intergenic} "
        "{params.promoter_args} "
        "{params.cpg_args} "
        "--log {log}"


if _te_active:

    # Final, transposon-subtracted genic/promoter/cpg_islands region BEDs
    # (what the existing boxplot pipeline actually consumes) -- trims only
    # the TE-overlapping bases out of each interval via `bedtools subtract`
    # (NOT `intersect -v`, which would drop a whole gene body just because
    # one small intronic TE fragment overlaps it). Re-sorted afterwards
    # since downstream `bedtools intersect -sorted` requires it.
    rule subtract_te_from_genic:
        input:
            raw=_region_outputs["genic"],
            te=resources.repeat_mask,
            chrom_sizes="resources/chrom_sizes.txt",
        output:
            "bed/genic.bed",
        log:
            "logs/resources/subtract_te_genic.log",
        conda:
            "../envs/deeptools.yaml"
        threads: 1
        resources:
            runtime=15,
            mem_mb=4000,
        shell:
            "bedtools subtract -a {input.raw} -b {input.te} | "
            "bedtools sort -i - -g {input.chrom_sizes} > {output} 2> {log}"

    # Intergenic = genome minus genes (generate_regions.py); transposons
    # are not genic-overlapping by definition here, so trim them out too.
    use rule subtract_te_from_genic as subtract_te_from_intergenic with:
        input:
            raw=_region_outputs["intergenic"],
            te=resources.repeat_mask,
            chrom_sizes="resources/chrom_sizes.txt",
        output:
            "bed/intergenic.bed",
        log:
            "logs/resources/subtract_te_intergenic.log",

    if resources.regulatory_gtf:

        use rule subtract_te_from_genic as subtract_te_from_promoter with:
            input:
                raw=_region_outputs["promoter"],
                te=resources.repeat_mask,
                chrom_sizes="resources/chrom_sizes.txt",
            output:
                "bed/promoter.bed",
            log:
                "logs/resources/subtract_te_promoter.log",

    if resources.cpg_islands:

        use rule subtract_te_from_genic as subtract_te_from_cpg_islands with:
            input:
                raw=_region_outputs["cpg_islands"],
                te=resources.repeat_mask,
                chrom_sizes="resources/chrom_sizes.txt",
            output:
                "bed/cpg_islands.bed",
            log:
                "logs/resources/subtract_te_cpg_islands.log",

    # Step 2 of the TE feature: drop any repeat element that overlaps a
    # gene body entirely (not a partial trim -- see subtract_te_from_genic
    # above for that), using the RAW (pre-subtraction) genic span so this
    # doesn't circularly depend on genic.bed already having TEs removed.
    # Shared by every configured TE class block below.
    rule filter_repeat_mask_nongenic:
        input:
            repeat_mask=resources.repeat_mask,
            genic=_region_outputs["genic"],
        output:
            "resources/repeat_mask_nongenic.bed",
        log:
            "logs/resources/filter_repeat_mask_nongenic.log",
        conda:
            "../envs/deeptools.yaml"
        threads: 1
        resources:
            runtime=30,
            mem_mb=8000,
        shell:
            "bedtools intersect -v -a {input.repeat_mask} -b {input.genic} "
            "> {output} 2> {log}"

    # Steps 3-4 of the TE feature, one rule per configured TE class block
    # (see te_class_blocks() in common.smk): subset to that class's
    # elements of at least the configured length, then optionally split
    # into families and, within those, subfamilies (repName prefix,
    # "_"-boundary matched) -- see generate_te_regions.py. Collisions
    # between class/family/subfamily names (within or across blocks, or
    # against standard/custom regions) are already forbidden by
    # validate_te_config() at Snakefile-parse time, so each rule's outputs
    # below are guaranteed globally unique.
    for _class_name, _block in te_class_blocks():
        _family_list = te_family_list(_block)
        _subfamilies = list(_block.get("subfamilies", []))
        _family_beds = {f: f"bed/{f}.bed" for f in _family_list}
        _subfamily_beds = {s: f"bed/{s}.bed" for s in _subfamilies}
        _family_args = " ".join(
            f"--family-bed {f}={p}" for f, p in _family_beds.items()
        )
        _subfamily_args = " ".join(
            f"--subfamily-bed {s}={p}" for s, p in _subfamily_beds.items()
        )

        rule:
            name:
                f"generate_te_regions_{_class_name}"
            input:
                repeat_mask="resources/repeat_mask_nongenic.bed",
            output:
                class_bed=f"bed/{_class_name}.bed",
                family=list(_family_beds.values()),
                subfamily=list(_subfamily_beds.values()),
            params:
                script=os.path.join(WORKFLOW_SCRIPTS, "generate_te_regions.py"),
                class_name=_class_name,
                min_length=_block["min_length"],
                family_args=_family_args,
                subfamily_args=_subfamily_args,
            log:
                f"logs/resources/generate_te_regions_{_class_name}.log",
            conda:
                "../envs/deeptools.yaml"
            threads: 1
            resources:
                runtime=30,
                mem_mb=8000,
            shell:
                "python {params.script} "
                "--repeat-mask {input.repeat_mask} "
                "--class-name {params.class_name} "
                "--min-length {params.min_length} "
                "--class-bed {output.class_bed} "
                "{params.family_args} "
                "{params.subfamily_args} "
                "--log {log}"


rule bismark_genome_preparation:
    input:
        fasta="resources/combined_genome.fa",
    output:
        directory("resources/Bisulfite_Genome"),
    log:
        "logs/resources/bismark_genome_preparation.log",
    threads: 40  # make sure to assign half of this to bismark
    resources:
        runtime=360,
        mem_mb=60000,
    conda:
        "../envs/bismark.yaml"
    script:
        "../scripts/bismark_genome_preparation.py"


# Annotate CpGs in the genome
# -----------------------------------------------------
rule find_cpgs:
    input:
        resources.filtered_fasta,
    output:
        "resources/cpg_sites.bed",
    params:
        script=os.path.join(WORKFLOW_SCRIPTS, "find_cpgs.py"),
    log:
        "logs/resources/find_cpgs.log",
    threads: 10
    resources:
        runtime=120,
        mem_mb=10000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "python {params.script} {input} {output} {log}"


# Create CpG probe BED file
# -----------------------------------------------------
rule create_cpg_probes:
    input:
        "resources/cpg_sites.bed",
    output:
        "resources/cpg_probes.bed",
    params:
        script=os.path.join(WORKFLOW_SCRIPTS, "create_cpg_probes.py"),
        # Number of CpGs per probe
        n=config["boxplot"]["cpg_n"],
    log:
        "logs/resources/create_cpg_probes.log",
    threads: 10
    resources:
        runtime=240,
        mem_mb=10000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "python {params.script} {input} {output} {params.n} {log}"


# Sort CpG probe BED file using chrom sizes
# -----------------------------------------------------
rule sort_cpg_probes:
    input:
        probes="resources/cpg_probes.bed",
        cs="resources/chrom_sizes.txt",
    output:
        "resources/cpg_probes_sorted.bed",
    log:
        "logs/resources/sort_cpg_probes.log",
    threads: 1
    resources:
        runtime=15,
        mem_mb=2000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bedtools sort -i {input.probes} -g {input.cs} > {output} 2> {log}"
