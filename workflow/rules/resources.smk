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


use rule get_genome_fasta as get_gtf with:
    output:
        resources.gtf,
    params:
        url=resources.gtf_url,
    log:
        "logs/resources/get_gtf.log",


use rule get_genome_fasta as get_regulatory_gtf with:
    output:
        resources.regulatory_gtf,
    params:
        url=resources.regulatory_gtf_url,
    log:
        "logs/resources/get_regulatory_gtf.log",


# Downloads the UCSC cpgIslandExt table (goldenPath database dump -- the same
# table hgTables serves under Track: CpG Islands) and converts it to a plain
# BED file, stripping the "chr" prefix and keeping only standard chromosomes
# to match the Ensembl-style contigs used elsewhere in this workflow.
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


# Downloads the UCSC rmsk table (goldenPath database dump -- the same table
# hgTables serves under Track: RepeatMasker), converts it to a plain BED
# file (chrom, start, end, repName, ., strand, repClass, repFamily),
# stripping the "chr" prefix and keeping only standard chromosomes to match
# the Ensembl-style contigs used elsewhere in this workflow, and drops
# non-transposable-element repeat classes/families listed in
# workflow/resources/nonTE_repClasses.txt.
rule prepare_repeat_mask:
    input:
        nonte="workflow/resources/nonTE_repClasses.txt",
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
        "'{{sub(/^chr/, \"\", $6); if ($6 ~ /^([0-9]+|X|Y|MT)$/) print $6, $7, $8, $11, \".\", $10, $13, $12}}' | "
        "grep -v -f {input.nonte} "
        "> {output} 2>> {log}"


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
    log:
        "logs/resources/find_cpgs.log",
    threads: 10
    resources:
        runtime=120,
        mem_mb=10000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "python workflow/scripts/find_cpgs.py {input} {output} {log}"


# Create CpG probe BED file
# -----------------------------------------------------
rule create_cpg_probes:
    input:
        "resources/cpg_sites.bed",
    output:
        "resources/cpg_probes.bed",
    params:
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
        "python workflow/scripts/create_cpg_probes.py {input} {output} {params.n} {log}"


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


"""
# Extend bed regions by a number of bases on each side
# -----------------------------------------------------
rule extend_bed_regions:
    input:
        bed="bed/{meta_region}.bed",
        cs="resources/chrom_sizes.txt",
    output:
        extended_bed="bed/{meta_region}_extended.bed",
    params:
        extend_with=config["metaplot"]["extend"],
    log:
        "logs/resources/extend_bed_{meta_region}.log",
    threads: 1
    resources:
        runtime=15,
        mem_mb=2000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        "grep -v 'random' {input.bed} | "
        "bedtools slop -i stdin "
        "-g {input.cs} "
        "-b {params.extend_with} > {output.extended_bed} 2> {log}"
"""
