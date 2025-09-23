import pandas as pd
import glob
import os
from scripts.resources import Resources
from snakemake.utils import validate
from snakemake.logging import logger

resources = Resources(config["genome"], config["ensembl_genome_build"])

# validate sample sheet and config file
# validate(samples, schema="schemas/samples.schema.yaml")
# validate(config, schema="schemas/config.schema.yaml")


def paired_end():
    fastq = glob.glob("reads/*.fastq.gz")

    if len(fastq) == 0:
        raise ValueError("No FASTQ (*.fastq.gz) files found in 'reads/' directory.")

    paired_end = all(("_R1_" in f or "_R2_" in f) for f in fastq)

    if paired_end:
        logger.info("Paired-end reads detected.")
    else:
        logger.info("Single-end reads detected.")

    return paired_end


def import_samples(paired_end):

    if paired_end:
        fastq = glob.glob("reads/*_R1_001.fastq.gz")
        samples = [f.split("/")[-1].replace("_R1_001.fastq.gz", "") for f in fastq]
    else:
        fastq = glob.glob("reads/*.fastq.gz")
        samples = [f.split("/")[-1].replace(".fastq.gz", "") for f in fastq]

    return samples


def conditions(csv):
    return list(set(csv["condition"]))


def dedup_input(wildcards):
    if PAIRED_END:
        return {
            "bam": "results/bismark/{wildcards.sample}/{wildcards.sample}_R1_bismark_bt2_pe.bam".format(
                wildcards=wildcards
            )
        }
    else:
        return {
            "bam": "results/bismark/{wildcards.sample}/{wildcards.sample}_bismark_bt2.bam".format(
                wildcards=wildcards
            )
        }
