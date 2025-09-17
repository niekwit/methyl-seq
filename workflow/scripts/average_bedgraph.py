"""
This script averages signal values across multiple compressed bedGraph files for each experimental condition.
It is designed to be used within a Snakemake workflow, utilizing Snakemake's input, output, and wildcards variables.
For each condition, it:
- Identifies all relevant bedGraph files (compressed with gzip).
- Decompresses them on the fly.
- Uses `bedtools unionbedg` to merge the files by genomic coordinates.
- Uses `awk` to compute the average signal value across replicates for each region.
- Writes the averaged bedGraph to the specified output file.
Logging is set up to record input parameters, processing steps, and executed commands for debugging purposes.
"""

import logging
from snakemake.shell import shell

# Set up logging
log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(log)],
)

# Snakemake variables
bg_files = snakemake.input["bg"]
condition = snakemake.wildcards["condition"]
output_file = snakemake.output["bg"]

# Log variables for debugging
logging.info(f"Input bedGraph files: {bg_files}")
logging.info(f"Condition: {condition}")
logging.info(f"Output bedGraph files: {output_file}")

# Create average bedGraph file
logging.info(f"Processing condition: {condition}")

# Each BedGraph file is compressed, so we need to decompress them on the fly
bg_input = " ".join([f"<(zcat {f})" for f in bg_files])

# Use awk to create average signal column
# escape braces for Snakemake + awk
replicate_count = len(bg_files)
awk_sum = "+".join([f"${i}" for i in range(4, 4 + replicate_count)])
awk = (
    "awk 'BEGIN{{OFS=\"\\t\"}} "
    "{{print $1, $2, $3, (" + awk_sum + ")/" + str(replicate_count) + "}}'"
)

# Create sort subcommand to filter out duplicate regions
# and merge them by averaging the signal
merge_sort = "LC_ALL=C sort -k1,1 -k2,2n | bedtools merge -c 4 -o mean -i - "

# Create full command
command = f"bedtools unionbedg -i {bg_input} | {awk} | {merge_sort} > {output_file}"
logging.info(f"Running command: {command}")

# Execute command
shell(command)
