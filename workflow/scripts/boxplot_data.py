"""

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

# Get data files
cpg_probe_file = snakemake.input["probes"]
methylation_file = snakemake.input["meth"]

# Get region
region = snakemake.wildcards["region"]

