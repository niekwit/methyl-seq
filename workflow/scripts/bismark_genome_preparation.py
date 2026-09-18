# Use half of the assigned threads in the Bismark command:
# Top and bottom strands are indexed separately

import glob
import os
from snakemake.shell import shell

# Get current working dir
cwd = os.getcwd()

# bismark_genome_preparation globs *.fa (then *.fa.gz, *.fasta, *.fasta.gz)
# directly in the given folder (non-recursive) and treats each match as a
# separate genome file, checking only that each file's FIRST FastA header
# is unique across files -- it does not recheck headers it encounters
# later within the same multi-FastA file. If more than one *.fa file ever
# ends up directly under resources/ (e.g. a leftover/stale one), a shared
# first-chromosome name between them raises exactly the "chromosome name
# ... already exists" error this is diagnosing. Log what's actually there
# right before invoking it.
fa_files = sorted(glob.glob("resources/*.fa"))
print(f"resources/*.fa before bismark_genome_preparation: {fa_files}")
for f in fa_files:
    with open(f) as fh:
        first_line = fh.readline().rstrip("\n")
    print(f"  {f}: first header = {first_line!r}")

# bismark_genome_preparation rejects --parallel < 2; with few cores
# available (e.g. CI's --cores 2), snakemake.threads can be scaled down low
# enough that halving it would go below that floor.
bismark_threads = max(2, int(snakemake.threads / 2))
command = (
    f"bismark_genome_preparation --verbose --parallel {bismark_threads} resources/"
)
print(command)
shell(command)
