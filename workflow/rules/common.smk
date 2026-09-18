import pandas as pd
import glob
import os
import re
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
    condition="[^/]+",


def targets():
    targets = [
        "results/plots/PCA.pdf",
        # "results/multiqc/multiqc_bismark.html",
        expand("results/bigwig/{condition}.bw", condition=CONDITIONS),
        expand("results/bigwig/coverage/{condition}.bw", condition=CONDITIONS),
        "results/plots/methylation_conversion_rate.csv",
        "results/plots/methylation_conversion_rate.pdf",
        "results/multiqc/multiqc_report.html",
    ]

    if config["boxplot"]["plot"]:
        targets.append("results/plots/boxplots.pdf")

    if te_regions():
        targets.append("results/plots/te_boxplots.pdf")

    if config["DMR"]["run"]:
        targets.extend(
            [
                "results/dmrs/hypomethylated_DMRs_annotated.tab",
                "results/dmrs/hypermethylated_DMRs_annotated.tab",
                "results/plots/dmrs/DMR_distance_to_TSS.pdf",
                "results/plots/dmrs/DMR_genomic_distribution.pdf",
                "results/plots/dmrs/DMR_volcano.pdf",
                "results/dmrs/all_methylation_tiles.rds",
                "results/dmrs/differential_methylation_tiles.rds",
                "results/dmrs/significant_differential_methylation_tiles.rds",
            ]
        )

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


def samples_in_condition(csv, condition):
    """
    Sample names belonging to a condition, via the exact sample/condition
    mapping in samples.csv -- not a name-prefix heuristic. A prefix check
    (sample.startswith(condition)) silently pulls in unrelated samples
    whenever one sample's name is a prefix of another's, e.g. condition
    "WT" wrongly matching a sample named "WT_3" whose own condition is the
    separate "WT_3".
    """
    return csv.loc[csv["condition"] == condition, "sample"].tolist()


def validate_reference_condition(reference_condition, conditions):
    """
    Checks that config DMR:reference_condition names an actual condition
    from samples.csv, so a typo fails fast at DAG-build time instead of
    surfacing later as a confusing methylKit/R error (or silently
    comparing against nothing, since dmr.R matches it via string
    detection against sample names).
    """
    if reference_condition not in conditions:
        raise ValueError(
            f"DMR reference_condition '{reference_condition}' is not one of "
            f"the conditions in config/samples.csv: {sorted(conditions)}"
        )


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
    """
    Standard genomic regions generated automatically by the generate_regions
    rule (see resources.smk), plus any optional extra custom regions the
    user defines under config boxplot:regions (name -> BED file path).
    """
    standard = ["whole_genome", "genic", "exon", "intron", "intergenic"]
    if resources.regulatory_gtf:
        standard.append("promoter")
    if resources.cpg_islands:
        standard.append("cpg_islands")

    extra = config["boxplot"].get("regions", None) or {}

    clashes = set(extra.keys()) & (set(standard) | set(te_regions()))
    if clashes:
        raise ValueError(
            f"config boxplot:regions name(s) {sorted(clashes)} clash with "
            "automatically generated standard/TE region(s) of the same "
            "name -- please rename them"
        )

    return standard + list(extra.keys())


_BOXPLOT_RESERVED = {"plot", "cpg_n", "min_reads", "regions"}


def te_class_blocks():
    """
    Ordered (class_name, block) pairs for every non-reserved key directly
    under config boxplot -- each such key is a TE class block, keyed by a
    repeat_mask.bed repClass value (e.g. LINE, LTR, SINE). Order matches
    config file order (dict insertion order), which is also plotting
    order in the combined TE boxplot.
    """
    return [
        (name, block)
        for name, block in config["boxplot"].items()
        if name not in _BOXPLOT_RESERVED
    ]


def te_family_list(block):
    """
    Parses a TE class block's `family` (an optional comma-separated
    string, e.g. "L1,L2") into an ordered, de-duplicated list. [] if
    `family` is absent.
    """
    raw = block.get("family")
    if not raw:
        return []
    seen = []
    for f in raw.split(","):
        f = f.strip()
        if f and f not in seen:
            seen.append(f)
    return seen


def te_regions():
    """
    Ordered list of every TE region name (class total, then family
    totals, then subfamily totals, per block) across all configured TE
    class blocks under config boxplot, in config order. A class block's
    mere presence enables it -- there is no separate plot flag per block.
    Returns [] when no TE class blocks are configured -- REGIONS,
    targets(), and the region wildcard_constraint are then completely
    unaffected, driving the separate combined TE boxplot figure
    (results/plots/te_boxplots.pdf) instead of the main one.
    """
    names = []
    for class_name, block in te_class_blocks():
        names.append(class_name)
        names.extend(te_family_list(block))
        names.extend(block.get("subfamilies", []))
    return names


def validate_te_config(te_regions, other_regions, resources):
    """
    Static (config-only) checks at Snakefile-parse time:
    1. TE class blocks require a genome with a RepeatMasker track
       (resources.repeat_mask is None e.g. for dm6).
    2. No TE region name (class/family/subfamily, across all blocks)
       collides with another TE region name, or with a standard/custom
       region name from regions().
    Whether a configured class/family/subfamily actually matches any
    elements in repeat_mask.bed can't be checked here (repeat_mask.bed is
    only materialised by a job at DAG-execution time) -- that hard-errors
    inside generate_te_regions.py instead, mirroring this function's
    relationship to validate_reference_condition() above.
    """
    if te_regions and resources.repeat_mask is None:
        raise ValueError(
            f"config boxplot has TE class block(s), but genome "
            f"'{resources.genome}' has no RepeatMasker track available "
            "for this pipeline (Resources.repeat_mask_url is not set) -- "
            "TE region analysis is not possible."
        )

    seen, dupes = set(), set()
    for name in te_regions:
        (dupes if name in seen else seen).add(name)
    if dupes:
        raise ValueError(
            f"TE region name(s) {sorted(dupes)} are produced by more than "
            "one boxplot class/family/subfamily entry -- region names must "
            "be unique across all TE class blocks; rename the colliding "
            "family/subfamily value(s) or class block key(s)."
        )

    clashes = set(te_regions) & set(other_regions)
    if clashes:
        raise ValueError(
            f"TE region name(s) {sorted(clashes)} clash with standard/"
            "custom boxplot region(s) of the same name -- please rename them."
        )


def meta_regions():
    regions = config["metaplot"].get("regions", None)
    # Get keys of regions dict
    if regions:
        return list(regions.keys())
    else:
        raise ValueError("No regions defined in config file under metaplot:regions")
