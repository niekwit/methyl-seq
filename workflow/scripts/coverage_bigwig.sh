#!/bin/bash
set -e # Exit on error

# Usage: script.sh input.cov.gz chrom.sizes temp_bedgraph output.bw

INPUT="$1"
SIZES="$2"
TEMP_BG="$3" # Snakemake will handle cleanup of this temporary file
OUTPUT="$4"

# Data structure of .cov.gz:
# <chrom> <start> <end> <methylated_percentage> <methylated_count> <unmethylated_count>

# Extract total depth ($5+$6), fix 0-length coordinates ($2, $2+1)
# Sort using LC_COLLATE=C for bedGraphToBigWig compatibility
zcat "$INPUT" | \
    awk 'OFS="\t" { print $1, $2, $2+1, $5+$6 }' | \
    LC_COLLATE=C sort -k1,1 -k2,2n > "$TEMP_BG"

# Convert to BigWig
bedGraphToBigWig "$TEMP_BG" "$SIZES" "$OUTPUT"