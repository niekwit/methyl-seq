# import basic packages
import pandas as pd
import glob
import os
from scripts.resources import Resources
from snakemake.utils import validate

resources = Resources(config["genome"], config["ensembl_genome_build"])

# Get current working dir
cwd = os.getcwd()

# validate sample sheet and config file
#validate(samples, schema="schemas/samples.schema.yaml")
#validate(config, schema="schemas/config.schema.yaml")

def import_samples():
    
    fastq = glob.glob("reads/*_R1_001.fastq.gz")
    samples = [f.split("/")[-1].replace("_R1_001.fastq.gz", "") for f in fastq]
    return samples