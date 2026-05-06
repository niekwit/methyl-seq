import pandas as pd
import glob
import os
from scripts.resources import Resources
from snakemake.utils import validate
from snakemake.logging import logger

resources = Resources(config["genome"], config["ensembl_genome_build"])


def check_cwd_no_spaces():
    cwd = os.path.abspath(os.getcwd())
    if " " in cwd:
        raise ValueError(
            f"Working directory path contains spaces: '{cwd}'\n"
            "Snakemake does not support paths with spaces. Please move the project to a path without spaces."
        )


# validate config file
validate(config, schema="../schemas/config.schema.yaml")


wildcard_constraints:
    # do not include slashes in sample and condition wildcards to avoid 
    # issues with rules that use these wildcards in file paths
    sample="[^/]+",
    condition="[^/]+"


def targets():
    targets = [
        "results/plots/PCA.pdf",
        # "results/multiqc/multiqc_bismark.html",
        expand("results/bigwig/{condition}.bw", condition=CONDITIONS),
        expand("results/bigwig/coverage/{condition}.bw", condition=CONDITIONS),
        "results/plots/methylation_conversion_rate.csv",
        "results/plots/methylation_conversion_rate.pdf",
        "results/preseq/library_complexity_summary.txt",
        "results/multiqc/multiqc_report.html",
    ]

    if config["boxplot"]["plot"]:
        targets.append("results/plots/boxplots.pdf")

    if config["DMR"]["run"]:
        targets.extend([
            "results/dmrs/hypomethylated_DMRs_annotated.tab",
            "results/dmrs/hypermethylated_DMRs_annotated.tab",
            "results/plots/dmrs/DMR_distance_to_TSS.pdf",
            "results/plots/dmrs/DMR_genomic_distribution.pdf",
            "results/plots/dmrs/DMR_volcano.pdf",
            "results/dmrs/all_methylation_tiles.rds",
            "results/dmrs/differential_methylation_tiles.rds",
            "results/dmrs/significant_differential_methylation_tiles.rds",
        ])

    return targets


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


def regions():
    regions = config["boxplot"].get("regions", None)
    # Get keys of regions dict
    if regions:
        # Create empty bed file for whole genome
        with open("bed/whole.genome.bed", "w") as f:
            f.write("")

        return ["whole.genome"] + list(regions.keys())
    else:
        raise ValueError("No regions defined in config file under boxplot:regions")


def meta_regions():
    regions = config["metaplot"].get("regions", None)
    # Get keys of regions dict
    if regions:
        return list(regions.keys())
    else:
        raise ValueError("No regions defined in config file under metaplot:regions")
