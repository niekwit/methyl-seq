#!/bin/bash

REGION=${snakemake_wildcards[region]}

if [ "$REGION" == "whole.genome" ]; then
	cp "${snakemake_input[probes]}" "${snakemake_output[region_probes]}"
else
	bedtools intersect -sorted -wa -a "${snakemake_input[probes]}" -b "${snakemake_input[regions]}" | \
	sort -k1,1 -k2,2n > "${snakemake_output[region_probes]}"
fi