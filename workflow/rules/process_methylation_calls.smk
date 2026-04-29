# Convert CpG methylation calls of OT to BED (merge replicates)
# -----------------------------------------------------
rule OT_methylation_calls_to_bed:
    input:
        cpgot=lambda wildcards: expand(
            "results/bismark/{sample}/CpG_OT_{sample}.deduplicated.txt.gz",
            sample=[s for s in SAMPLES if s.startswith(wildcards.condition)],
        ),
    output:
        bed=temp("results/bed/CpG_OT_{condition}.bed"),
    log:
        "logs/methylation_calls_to_bed/OT_{condition}.log",
    threads: 4
    resources:
        runtime=240,
    conda:
        "../envs/bismark.yaml"
    shell:
        "zcat {input.cpgot} | "
        "grep -vE 'phage_lambda|plasmid_puc19c' | "
        'awk -v OFS="\t" \'{{print $3, $4, $4+1, $1, $5, "+"}}\' > {output.bed}'


# Convert CpG methylation calls of OB to BED (merge replicates)
# -----------------------------------------------------
rule OB_methylation_calls_to_bed:
    input:
        cpgob=lambda wildcards: expand(
            "results/bismark/{sample}/CpG_OB_{sample}.deduplicated.txt.gz",
            sample=[s for s in SAMPLES if s.startswith(wildcards.condition)],
        ),
    output:
        bed=temp("results/bed/CpG_OB_{condition}.bed"),
    log:
        "logs/methylation_calls_to_bed/OB_{condition}.log",
    threads: 4
    resources:
        runtime=240,
    conda:
        "../envs/bismark.yaml"
    shell:
        "zcat {input.cpgob} | "
        "grep -vE 'phage_lambda|plasmid_puc19c' | "
        'awk -v OFS="\t" \'{{print $3, $4, $4+1, $1, $5, "-"}}\' > {output.bed}'


# Merge CpG methylation calls of OT and OB to a single BED file
# -----------------------------------------------------
"""
Performs a memory-efficient, genome-ordered sort and merge of BED files.

This command avoids the high RAM overhead of `bedtools sort` on large datasets 
(>20GB) by using a streaming Unix sort strategy:
1. Awk creates a numeric 'rank' lookup from the chrom_sizes file.
2. Input BED lines are prepended with this rank, effectively mapping 
   non-alphanumeric chromosome orders (e.g., MT, X, Y) to 
   integers.
3. Unix sort performs a parallelized, disk-buffered numeric sort on the 
   rank and start positions.
4. The temporary rank prefix is removed to restore standard BED format.

Constraints:
- Chromosomes in the BED file missing from the chrom_sizes file will 
  be filtered out.
- Memory usage is strictly controlled by the -S flag.
"""


rule merge_strand_methylation_calls_to_bed:
    input:
        ot="results/bed/CpG_OT_{condition}.bed",
        ob="results/bed/CpG_OB_{condition}.bed",
        cs="resources/chrom_sizes.txt",
    output:
        bed="results/bed/CpG_merged_{condition}.bed",
    log:
        "logs/methylation_calls_to_bed/merge_{condition}.log",
    params:
        tmpdir=config["temp_dir"],
    threads: 8
    resources:
        runtime=120,
        mem_mb=12000,
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        (cat {input.ot} {input.ob} | \
        awk 'NR==FNR {{ rank[$1]=NR; next }} ($1 in rank) {{ print rank[$1], $0 }}' {input.cs} - | \
        sort -k1,1n -k3,3n \
            --parallel={threads} \
            -S {resources.mem_mb}M \
            -T {params.tmpdir} | \
        cut -d' ' -f2- > {output.bed}) 2> {log}
        """
