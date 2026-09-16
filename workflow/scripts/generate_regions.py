#!/usr/bin/env python3

"""
Generates BED files of genomic regions from an Ensembl-style genome GTF and
(optionally) an Ensembl regulatory-build GFF3 file.

From the genome GTF:
    - whole genome : one interval per contig, spanning 0 to the highest
                      coordinate annotated on that contig in the GTF. This is
                      an approximation of chromosome length derived purely
                      from the GTF (no FASTA/.fai is read) -- use
                      resources/chrom_sizes.txt instead if exact assembly
                      lengths are required.
    - genic regions: one interval per feature where column 3 == "gene".

From the regulatory GTF/GFF3 (optional):
    - promoters    : one interval per feature where column 3 == "promoter".

An existing CpG island BED file (e.g. downloaded from UCSC) can also be
passed in via --cpg-island-bed; it is not derived from the GTF and is only
ever cleaned of blacklist overlaps (see below), never regenerated.

Any of the produced/passed-in BED files can optionally be cleaned of regions
that overlap a user-supplied "exclude" BED file (e.g. a blacklist), via
`bedtools intersect -v`.

Usage:
    python generate_regions.py \\
        --gtf genome.gtf.gz \\
        --whole-genome-bed whole_genome.bed \\
        --genic-bed genic.bed \\
        [--regulatory-gtf regulatory_features.gff3.gz --promoter-bed promoters.bed] \\
        [--cpg-island-bed cpg_islands.bed] \\
        [--exclude-bed blacklist.bed] \\
        --log generate_regions.log

Requirements:
    - bedtools on PATH (only needed when --exclude-bed is given)
"""

import argparse
import gzip
import logging
import re

import pybedtools


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def feature_name(feature, attributes):
    """Best-effort ID for a GTF/GFF3 feature, for the BED 'name' column."""
    # GTF style: key "value";
    for key in ("gene_id", "gene_name"):
        m = re.search(rf'{key} "([^"]+)"', attributes)
        if m:
            return m.group(1)
    # GFF3 style: key=value;
    m = re.search(r"ID=([^;]+)", attributes)
    if m:
        return m.group(1).rsplit(":", 1)[-1]
    return feature


def iter_gtf(path):
    """Yields (chrom, feature, start0, end, name) for every non-comment line."""
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, feature, start, end = fields[0:5]
            attributes = fields[8]
            # GTF/GFF3 coordinates are 1-based inclusive -> BED is 0-based half-open
            yield chrom, feature, int(start) - 1, int(end), feature_name(
                feature, attributes
            )


def whole_genome_bed(gtf, out_bed):
    logging.info(
        f"Deriving whole-genome intervals (max annotated coordinate per contig) from {gtf}"
    )
    max_end = {}
    for chrom, _feature, _start0, end, _name in iter_gtf(gtf):
        if end > max_end.get(chrom, 0):
            max_end[chrom] = end

    with open(out_bed, "wt") as fh:
        for chrom in sorted(max_end):
            fh.write(f"{chrom}\t0\t{max_end[chrom]}\n")

    logging.info(f"Wrote {len(max_end)} contigs to {out_bed}")


def feature_bed(gtf, feature, out_bed):
    logging.info(f"Extracting '{feature}' features from {gtf}")
    records = [
        (chrom, start0, end, name)
        for chrom, feat, start0, end, name in iter_gtf(gtf)
        if feat == feature
    ]
    records.sort(key=lambda r: (r[0], r[1]))

    with open(out_bed, "wt") as fh:
        for chrom, start0, end, name in records:
            fh.write(f"{chrom}\t{start0}\t{end}\t{name}\n")

    logging.info(f"Wrote {len(records)} '{feature}' intervals to {out_bed}")


def exclude_overlaps(bed_path, exclude_bed):
    logging.info(f"Removing regions in {bed_path} that overlap {exclude_bed}")
    n_before = sum(1 for _ in open(bed_path))

    kept = pybedtools.BedTool(bed_path).intersect(exclude_bed, v=True)
    kept.moveto(bed_path)

    n_after = sum(1 for _ in open(bed_path))
    logging.info(f"{bed_path}: {n_before} -> {n_after} intervals after exclusion")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gtf", required=True, help="Genome GTF file (.gtf or .gtf.gz)"
    )
    parser.add_argument(
        "--whole-genome-bed",
        required=True,
        help="Output BED for whole-genome intervals",
    )
    parser.add_argument(
        "--genic-bed",
        required=True,
        help="Output BED for genic ('gene' feature) regions",
    )
    parser.add_argument(
        "--regulatory-gtf",
        help="Regulatory build GTF/GFF3 file (.gff3 or .gff3.gz), for promoter regions",
    )
    parser.add_argument(
        "--promoter-bed", help="Output BED for 'promoter' feature regions"
    )
    parser.add_argument(
        "--cpg-island-bed",
        help=(
            "Existing CpG island BED file (e.g. from UCSC). Not regenerated -- "
            "only cleaned of --exclude-bed overlaps, in place, when given"
        ),
    )
    parser.add_argument(
        "--exclude-bed",
        help="BED file of regions to exclude from every output/passed-in BED (via bedtools intersect -v)",
    )
    parser.add_argument("--log", help="Log file (default: stderr)")
    args = parser.parse_args()

    if bool(args.regulatory_gtf) != bool(args.promoter_bed):
        parser.error("--regulatory-gtf and --promoter-bed must be given together")

    return args


def main():
    args = parse_args()

    logging.basicConfig(
        format="%(levelname)s:%(asctime)s: %(message)s",
        level=logging.DEBUG,
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=(
            [logging.FileHandler(args.log)] if args.log else [logging.StreamHandler()]
        ),
    )

    whole_genome_bed(args.gtf, args.whole_genome_bed)
    feature_bed(args.gtf, "gene", args.genic_bed)

    out_beds = [args.whole_genome_bed, args.genic_bed]

    if args.regulatory_gtf:
        feature_bed(args.regulatory_gtf, "promoter", args.promoter_bed)
        out_beds.append(args.promoter_bed)

    if args.cpg_island_bed:
        out_beds.append(args.cpg_island_bed)

    if args.exclude_bed:
        for bed in out_beds:
            exclude_overlaps(bed, args.exclude_bed)

    logging.info("Done.")


if __name__ == "__main__":
    main()
