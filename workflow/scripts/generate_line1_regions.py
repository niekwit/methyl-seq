#!/usr/bin/env python3

"""
Generates LINE1 (L1 family) region BED files from a gene-body-excluded
RepeatMasker BED (see the filter_repeat_mask_nongenic rule in
resources.smk, which produces this script's --repeat-mask input via
`bedtools intersect -v` against the raw genic region).

From --repeat-mask (columns: chrom, start, end, repName, score, strand,
repClass, repFamily -- see prepare_repeat_mask in resources.smk):
    - LINE1 : every element with repFamily == "L1" (repClass == "LINE"),
              length >= --min-length.
    - subfamily subsets : for each --subfamily-bed NAME=PATH, the subset of
              the above whose repName matches NAME exactly or NAME followed
              by "_" (e.g. "L1MdA" matches "L1MdA"/"L1MdA_I"/"L1MdA_II"/...,
              but not "L1MdFanc_I" when NAME is "L1MdF"). A subfamily that
              matches zero elements is a hard error, since that almost
              always means a typo or an overly strict --min-length.

Usage:
    python generate_line1_regions.py \\
        --repeat-mask repeat_mask_nongenic.bed \\
        --min-length 5000 \\
        --line1-bed LINE1.bed \\
        --subfamily-bed L1MdA=L1MdA.bed \\
        --subfamily-bed L1MdTf=L1MdTf.bed \\
        --log generate_line1_regions.log
"""

import argparse
import logging
import re

import pandas as pd

REPEAT_MASK_COLUMNS = [
    "chrom",
    "start",
    "end",
    "repName",
    "score",
    "strand",
    "repClass",
    "repFamily",
]


def load_line1(repeat_mask, min_length):
    logging.info(f"Loading {repeat_mask}")
    df = pd.read_csv(
        repeat_mask,
        sep="\t",
        header=None,
        names=REPEAT_MASK_COLUMNS,
        dtype={"chrom": str},
    )
    logging.info(f"{len(df)} gene-body-excluded repeat elements total")

    line1 = df[(df["repFamily"] == "L1") & (df["repClass"] == "LINE")].copy()
    logging.info(f"{len(line1)} LINE/L1 elements before length filtering")

    line1["length"] = line1["end"] - line1["start"]
    line1 = line1[line1["length"] >= min_length]
    line1 = line1.sort_values(["chrom", "start"])
    logging.info(f"{len(line1)} LINE/L1 elements with length >= {min_length}")

    return line1


def write_bed4(df, out_bed, name_col="repName"):
    df[["chrom", "start", "end", name_col]].to_csv(
        out_bed, sep="\t", header=False, index=False
    )


def subfamily_subset(line1, name):
    """repName == name, or name followed by "_" -- not a raw substring match."""
    pattern = re.compile(rf"^{re.escape(name)}(_.*)?$")
    return line1[line1["repName"].str.match(pattern)]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repeat-mask",
        required=True,
        help="Gene-body-excluded RepeatMasker BED (8 columns, see prepare_repeat_mask)",
    )
    parser.add_argument(
        "--min-length",
        required=True,
        type=int,
        help="Minimum LINE1 element length (bp), applied after gene-body exclusion",
    )
    parser.add_argument(
        "--line1-bed", required=True, help="Output BED for all filtered LINE1 elements"
    )
    parser.add_argument(
        "--subfamily-bed",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help=(
            "LINE1 subfamily to subset and its output BED path, as NAME=PATH. "
            "Repeatable, one per configured subfamily, in plotting order."
        ),
    )
    parser.add_argument("--log", help="Log file (default: stderr)")
    args = parser.parse_args()

    subfamily_beds = []
    for value in args.subfamily_bed:
        if "=" not in value:
            parser.error(f"--subfamily-bed must be NAME=PATH, got: {value!r}")
        name, path = value.split("=", 1)
        subfamily_beds.append((name, path))
    args.subfamily_beds = subfamily_beds

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

    line1 = load_line1(args.repeat_mask, args.min_length)
    write_bed4(line1, args.line1_bed)
    logging.info(f"Wrote {len(line1)} LINE1 intervals to {args.line1_bed}")

    for name, out_bed in args.subfamily_beds:
        subset = subfamily_subset(line1, name)
        if subset.empty:
            raise SystemExit(
                f"ERROR: configured LINE1 subfamily '{name}' matched 0 elements "
                f"in {args.repeat_mask} after gene-body exclusion and "
                f"min_length={args.min_length} filtering. Check '{name}' against "
                "repName values there, or lower boxplot:LINE1:min_length."
            )
        write_bed4(subset, out_bed)
        logging.info(f"Wrote {len(subset)} '{name}' intervals to {out_bed}")

    logging.info("Done.")


if __name__ == "__main__":
    main()
