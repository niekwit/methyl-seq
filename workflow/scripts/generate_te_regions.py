#!/usr/bin/env python3

"""
Generates region BED files for one transposable-element (TE) class from a
gene-body-excluded RepeatMasker BED (see the filter_repeat_mask_nongenic
rule in resources.smk, which produces this script's --repeat-mask input
via `bedtools intersect -v` against the raw genic region).

From --repeat-mask (columns: chrom, start, end, repName, score, strand,
repClass, repFamily -- see prepare_repeat_mask in resources.smk):
    - class total  : every element with repClass == --class-name, length
                      >= --min-length.
    - family subsets  : for each --family-bed NAME=PATH, the subset of the
                      class total whose repFamily == NAME.
    - subfamily subsets : for each --subfamily-bed NAME=PATH, the subset of
                      the class total whose repName matches NAME exactly
                      or NAME followed by "_" (e.g. "L1MdA" matches
                      "L1MdA"/"L1MdA_I"/"L1MdA_II"/..., but not
                      "L1MdFanc_I" when NAME is "L1MdF"), AND whose
                      repFamily is one of the names given via
                      --family-bed in this same invocation.

A class, family, or subfamily filter that matches zero elements is a hard
error, since that almost always means a typo, a case mismatch (repClass/
repFamily/repName values are matched case-sensitively), or an overly
strict --min-length.

Usage:
    python generate_te_regions.py \\
        --repeat-mask repeat_mask_nongenic.bed \\
        --class-name LINE \\
        --min-length 6000 \\
        --class-bed LINE.bed \\
        --family-bed L1=L1.bed \\
        --family-bed L2=L2.bed \\
        --subfamily-bed L1MdA=L1MdA.bed \\
        --subfamily-bed L1MdF=L1MdF.bed \\
        --log generate_te_regions.log
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


def load_class(repeat_mask, class_name, min_length):
    logging.info(f"Loading {repeat_mask}")
    df = pd.read_csv(
        repeat_mask,
        sep="\t",
        header=None,
        names=REPEAT_MASK_COLUMNS,
        dtype={"chrom": str},
    )
    logging.info(f"{len(df)} gene-body-excluded repeat elements total")

    cls = df[df["repClass"] == class_name].copy()
    if cls.empty:
        raise SystemExit(
            f"ERROR: TE class '{class_name}' matched 0 elements in "
            f"{repeat_mask}. Check the name against repeat_mask.bed's "
            "repClass values (case-sensitive)."
        )
    logging.info(f"{len(cls)} '{class_name}' elements before length filtering")

    cls["length"] = cls["end"] - cls["start"]
    cls = cls[cls["length"] >= min_length]
    if cls.empty:
        raise SystemExit(
            f"ERROR: TE class '{class_name}' matched 0 elements with "
            f"length >= {min_length} in {repeat_mask}. Lower "
            f"boxplot:{class_name}:min_length, or check the value."
        )
    cls = cls.sort_values(["chrom", "start"])
    logging.info(f"{len(cls)} '{class_name}' elements with length >= {min_length}")

    return cls


def write_bed4(df, out_bed, name_col="repName"):
    df[["chrom", "start", "end", name_col]].to_csv(
        out_bed, sep="\t", header=False, index=False
    )


def family_subset(cls, name, class_name, min_length, src):
    subset = cls[cls["repFamily"] == name]
    if subset.empty:
        raise SystemExit(
            f"ERROR: configured family '{name}' (class '{class_name}') "
            f"matched 0 elements in {src} after gene-body exclusion and "
            f"min_length={min_length} filtering. Check '{name}' against "
            "repFamily values there (case-sensitive)."
        )
    return subset


def subfamily_subset(cls, name, family_names, class_name, min_length, src):
    """repName == name, or name followed by "_" -- not a raw substring
    match -- restricted to repFamily in family_names."""
    pattern = re.compile(rf"^{re.escape(name)}(_.*)?$")
    subset = cls[
        cls["repName"].str.match(pattern) & cls["repFamily"].isin(family_names)
    ]
    if subset.empty:
        raise SystemExit(
            f"ERROR: configured subfamily '{name}' (class '{class_name}', "
            f"family filter {sorted(family_names)}) matched 0 elements in "
            f"{src} after gene-body exclusion and min_length={min_length} "
            "filtering. Check it against repName values there, or whether "
            "it actually belongs to one of the configured families."
        )
    return subset


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repeat-mask",
        required=True,
        help="Gene-body-excluded RepeatMasker BED (8 columns, see prepare_repeat_mask)",
    )
    parser.add_argument(
        "--class-name",
        required=True,
        help="TE class to select (repeat_mask.bed repClass value, e.g. LINE)",
    )
    parser.add_argument(
        "--min-length",
        required=True,
        type=int,
        help="Minimum element length (bp), applied after gene-body exclusion",
    )
    parser.add_argument(
        "--class-bed", required=True, help="Output BED for all filtered class elements"
    )
    parser.add_argument(
        "--family-bed",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help=(
            "Family to subset and its output BED path, as NAME=PATH. "
            "Repeatable, one per configured family, in plotting order."
        ),
    )
    parser.add_argument(
        "--subfamily-bed",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help=(
            "Subfamily to subset and its output BED path, as NAME=PATH. "
            "Repeatable, one per configured subfamily, in plotting order. "
            "Requires at least one --family-bed."
        ),
    )
    parser.add_argument("--log", help="Log file (default: stderr)")
    args = parser.parse_args()

    def parse_name_path_pairs(flag, values):
        pairs = []
        for value in values:
            if "=" not in value:
                parser.error(f"{flag} must be NAME=PATH, got: {value!r}")
            name, path = value.split("=", 1)
            pairs.append((name, path))
        return pairs

    args.family_beds = parse_name_path_pairs("--family-bed", args.family_bed)
    args.subfamily_beds = parse_name_path_pairs("--subfamily-bed", args.subfamily_bed)

    if args.subfamily_beds and not args.family_beds:
        parser.error("--subfamily-bed requires at least one --family-bed")

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

    cls = load_class(args.repeat_mask, args.class_name, args.min_length)
    write_bed4(cls, args.class_bed)
    logging.info(f"Wrote {len(cls)} '{args.class_name}' intervals to {args.class_bed}")

    for name, out_bed in args.family_beds:
        subset = family_subset(
            cls, name, args.class_name, args.min_length, args.repeat_mask
        )
        write_bed4(subset, out_bed)
        logging.info(f"Wrote {len(subset)} '{name}' intervals to {out_bed}")

    family_names = {name for name, _ in args.family_beds}
    for name, out_bed in args.subfamily_beds:
        subset = subfamily_subset(
            cls, name, family_names, args.class_name, args.min_length, args.repeat_mask
        )
        write_bed4(subset, out_bed)
        logging.info(f"Wrote {len(subset)} '{name}' intervals to {out_bed}")

    logging.info("Done.")


if __name__ == "__main__":
    main()
