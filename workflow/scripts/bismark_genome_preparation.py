# Use half of the assigned threads in the Bismark command:
# Top and bottom strands are indexed separately

import os
from snakemake.shell import shell

# Get current working dir
cwd = os.getcwd()

bismark_threads = int(snakemake.threads / 2)
command = f"bismark_genome_preparation --verbose --parallel {bismark_threads} resources/"
print(command)
shell(command)