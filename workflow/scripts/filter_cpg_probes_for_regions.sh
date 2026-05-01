#!/bin/bash

REGION=${snakemake_wildcards[region]}

if [ "$REGION" == "whole.genome" ]; then
	# Filter to chromosomes present in chrom_sizes and sort to match genome order,
	# so bedtools intersect -sorted -g is satisfied downstream
	grep -Fwf <(cut -f1 "${snakemake_input[chrom_sizes]}") "${snakemake_input[probes]}" | \
		bedtools sort -g "${snakemake_input[chrom_sizes]}" > "${snakemake_output[region_probes]}"
else
	bedtools intersect -sorted -wa \
		-a "${snakemake_input[probes]}" \
		-b "${snakemake_input[regions]}" \
		-g "${snakemake_input[chrom_sizes]}" | \
		sort -k1,1 -k2,2n > "${snakemake_output[region_probes]}"
fi