#!/usr/bin/env python3

"""
Generates BED files of genomic regions from a chrom_sizes file, an
Ensembl-style genome GTF, and (optionally) an Ensembl regulatory-build GFF3
file.

From --chrom-sizes (e.g. resources/chrom_sizes.txt, cut from a FASTA .fai):
    - whole genome : one interval per contig, spanning 0 to its exact length.

From the genome GTF:
    - genic regions: one interval per feature where column 3 == "gene".
    - exon regions : unique intervals from column 3 == "exon", de-duplicated
                      by coordinates (transcripts of the same gene commonly
                      share exons).
    - intron regions: one interval per gap between consecutive exons of the
                       same transcript, derived via a per-transcript
                       groupby/shift over all exon rows (vectorised, no
                       per-transcript Python loop -- see intron_bed()), also
                       de-duplicated by coordinates.
    - intergenic regions: whole-genome intervals minus genic regions minus
                       --exclude-bed (if given), via `bedtools subtract`.

From the regulatory GTF/GFF3 (optional):
    - promoters    : one interval per feature where column 3 == "promoter".

An existing CpG island BED file (e.g. downloaded from UCSC) can also be
copied through to --cpg-island-bed-out via --cpg-island-bed; it is not
derived from the GTF and is only ever copied and (optionally) cleaned of
blacklist overlaps (see below), never regenerated.

Any of the produced/passed-in BED files can optionally be cleaned of regions
that overlap a user-supplied "exclude" BED file (e.g. a blacklist), via
`bedtools intersect -v`.

Usage:
    python generate_regions.py \\
        --gtf genome.gtf.gz \\
        --chrom-sizes chrom_sizes.txt \\
        --whole-genome-bed whole_genome.bed \\
        --genic-bed genic.bed \\
        --exon-bed exon.bed \\
        --intron-bed intron.bed \\
        --intergenic-bed intergenic.bed \\
        [--regulatory-gtf regulatory_features.gff3.gz --promoter-bed promoters.bed] \\
        [--cpg-island-bed cpg_islands.bed --cpg-island-bed-out cpg_islands_out.bed] \\
        [--exclude-bed blacklist.bed] \\
        --log generate_regions.log

Requirements:
    - bedtools on PATH (only needed when --exclude-bed is given)
"""

import argparse
import gzip
import logging
import re
import shutil

import pandas as pd
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


def exon_dataframe(gtf):
    """Loads exon rows (chrom, start0, end, transcript_id) from a GTF into a DataFrame."""
    rows = []
    with open_maybe_gzip(gtf) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, start, end, attributes = fields[0], fields[3], fields[4], fields[8]
            m = re.search(r'transcript_id "([^"]+)"', attributes)
            if not m:
                continue
            # GTF coordinates are 1-based inclusive -> BED is 0-based half-open
            rows.append((chrom, int(start) - 1, int(end), m.group(1)))

    return pd.DataFrame(rows, columns=["chrom", "start0", "end", "transcript_id"])


def exon_bed(gtf, out_bed):
    """Writes unique exon intervals (de-duplicated by coordinates) from a GTF."""
    logging.info(f"Extracting unique exon regions from {gtf}")
    exons = exon_dataframe(gtf)
    exons = exons.drop_duplicates(subset=["chrom", "start0", "end"])
    exons = exons.sort_values(["chrom", "start0"])

    exons.to_csv(out_bed, sep="\t", header=False, index=False)

    logging.info(f"Wrote {len(exons)} unique exon intervals to {out_bed}")


def intron_bed(gtf, out_bed):
    """
    Derive intron intervals from exon features via a per-transcript
    groupby/shift, instead of looping over transcripts in Python: the loop
    re-scans the exon table per transcript (O(transcripts x exons)), which
    does not scale to a full genome GTF.
    """
    logging.info(f"Deriving intron regions from exons in {gtf}")
    exons = exon_dataframe(gtf).sort_values(["transcript_id", "start0"])

    grouped = exons.groupby("transcript_id", sort=False)
    next_start = grouped["start0"].shift(-1)
    has_next = next_start.notna()

    introns = pd.DataFrame(
        {
            "chrom": exons.loc[has_next, "chrom"],
            "start0": exons.loc[has_next, "end"],
            "end": next_start[has_next].astype(exons["start0"].dtype),
            "name": exons.loc[has_next, "transcript_id"],
        }
    )
    # Drop non-positive-length gaps from overlapping/adjacent exon annotations
    introns = introns[introns["end"] > introns["start0"]]
    # Transcripts of the same gene (or overlapping genes) commonly share
    # identical intron coordinates; collapse those to one interval instead
    # of emitting one row per transcript that has it
    introns = introns.drop_duplicates(subset=["chrom", "start0", "end"])
    introns = introns.sort_values(["chrom", "start0"])

    introns.to_csv(out_bed, sep="\t", header=False, index=False)

    logging.info(f"Wrote {len(introns)} intron intervals to {out_bed}")


def whole_genome_bed(chrom_sizes, out_bed):
    logging.info(f"Deriving whole-genome intervals from {chrom_sizes}")
    n = 0
    with open(chrom_sizes) as fh, open(out_bed, "wt") as out:
        for line in fh:
            if not line.strip():
                continue
            chrom, size = line.rstrip("\n").split("\t")[:2]
            out.write(f"{chrom}\t0\t{size}\n")
            n += 1

    logging.info(f"Wrote {n} contigs to {out_bed}")


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


def intergenic_bed(whole_genome_bed, genic_bed, exclude_bed, out_bed):
    """Whole-genome intervals minus genic regions minus exclude_bed (if given)."""
    logging.info(
        f"Deriving intergenic regions from {whole_genome_bed} minus {genic_bed}"
        + (f" minus {exclude_bed}" if exclude_bed else "")
    )
    regions = pybedtools.BedTool(whole_genome_bed).subtract(genic_bed)
    if exclude_bed:
        regions = regions.subtract(exclude_bed)
    regions.moveto(out_bed)

    n = sum(1 for _ in open(out_bed))
    logging.info(f"Wrote {n} intergenic intervals to {out_bed}")


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
        "--chrom-sizes",
        required=True,
        help="Two-column chrom/size file (e.g. resources/chrom_sizes.txt), for whole-genome intervals",
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
        "--exon-bed",
        required=True,
        help="Output BED for unique exon regions",
    )
    parser.add_argument(
        "--intron-bed",
        required=True,
        help="Output BED for intron regions (gaps between consecutive exons of a transcript)",
    )
    parser.add_argument(
        "--intergenic-bed",
        required=True,
        help="Output BED for intergenic regions (whole genome minus genic regions minus --exclude-bed)",
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
            "Existing CpG island BED file (e.g. from UCSC), to copy to "
            "--cpg-island-bed-out. Not regenerated -- only copied and "
            "(optionally) cleaned of --exclude-bed overlaps"
        ),
    )
    parser.add_argument(
        "--cpg-island-bed-out",
        help="Output path for --cpg-island-bed (required together with --cpg-island-bed)",
    )
    parser.add_argument(
        "--exclude-bed",
        help="BED file of regions to exclude from every output/passed-in BED (via bedtools intersect -v)",
    )
    parser.add_argument("--log", help="Log file (default: stderr)")
    args = parser.parse_args()

    if bool(args.regulatory_gtf) != bool(args.promoter_bed):
        parser.error("--regulatory-gtf and --promoter-bed must be given together")

    if bool(args.cpg_island_bed) != bool(args.cpg_island_bed_out):
        parser.error("--cpg-island-bed and --cpg-island-bed-out must be given together")

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

    whole_genome_bed(args.chrom_sizes, args.whole_genome_bed)
    feature_bed(args.gtf, "gene", args.genic_bed)
    exon_bed(args.gtf, args.exon_bed)
    intron_bed(args.gtf, args.intron_bed)
    # Built directly from the not-yet-cleaned whole-genome/genic BEDs below,
    # with its own --exclude-bed subtraction -- not added to out_beds, since
    # it is already blacklist-clean by construction
    intergenic_bed(
        args.whole_genome_bed, args.genic_bed, args.exclude_bed, args.intergenic_bed
    )

    out_beds = [
        args.whole_genome_bed,
        args.genic_bed,
        args.exon_bed,
        args.intron_bed,
    ]

    if args.regulatory_gtf:
        feature_bed(args.regulatory_gtf, "promoter", args.promoter_bed)
        out_beds.append(args.promoter_bed)

    if args.cpg_island_bed:
        logging.info(f"Copying {args.cpg_island_bed} to {args.cpg_island_bed_out}")
        shutil.copyfile(args.cpg_island_bed, args.cpg_island_bed_out)
        out_beds.append(args.cpg_island_bed_out)

    if args.exclude_bed:
        for bed in out_beds:
            exclude_overlaps(bed, args.exclude_bed)

    logging.info("Done.")


if __name__ == "__main__":
    main()
