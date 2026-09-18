# Use half of the assigned threads in the Bismark command:
# Top and bottom strands are indexed separately

import os
from snakemake.shell import shell

# Get current working dir
cwd = os.getcwd()

# bismark_genome_preparation rejects --parallel < 2; with few cores
# available (e.g. CI's --cores 2), snakemake.threads can be scaled down low
# enough that halving it would go below that floor.
bismark_threads = max(2, int(snakemake.threads / 2))
command = (
    f"bismark_genome_preparation --verbose --parallel {bismark_threads} resources/"
)
print(command)
shell(command)
