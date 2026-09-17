#!/usr/bin/env python3

"""
Build a tiny methyl-seq test dataset for CI from real Bismark output.

It keeps every read (both mates) that Bismark aligned to a single small
locus (+ flanks), plus a capped number of spike-in control reads, and writes
them back out as FASTQ under .test/reads/. Optionally it also carves the
matching mini genome FASTA and mini annotation tracks (GTF, regulatory
GFF3, CpG islands, RepeatMasker) out of the full reference, so the "test"
genome branch in resources.py has something to serve for every code path.

Nothing is subsampled by default: coverage depth over the locus is preserved,
only the genomic breadth is reduced, so methylKit/boxplots stay meaningful.

BAM reading/scanning is done with pysam rather than shelling out to
`samtools view`/`awk`. These Bismark dedup BAMs aren't coordinate-indexed
(no random-access shortcut), so this is still a full linear pass either
way, and pysam's per-record Python object overhead does make it slower
than `samtools view | awk` in C -- roughly 1.5-2x on real ~5GB samples,
even after the two optimizations already applied here (multi-threaded
BGZF decompression, and matching control contigs once up front instead
of a regex call per read). That's an accepted, deliberate tradeoff: this
script only runs manually and rarely (regenerating committed test
fixtures), never in CI, so the extra wall-clock time doesn't matter in
practice and isn't worth splitting the codebase across two languages for.
The final FASTQ reconstruction (`samtools collate` + `samtools fastq`) is
left to samtools itself regardless: getting Bismark's paired-mate
ordering and strand bookkeeping exactly right is what that tool is for,
and hand-rolling it would risk silently corrupting the checked-in test
fixtures for no real benefit.

Requirements: pysam, samtools >= 1.12 (needs `view -N` / --qname-file).
"""

import argparse
import glob
import gzip
import os
import re
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import pysam

OUT_DIR = Path(__file__).resolve().parent
READS_DIR = OUT_DIR / "reads"
RES_DIR = OUT_DIR / "resources"

# Target region (mm39). The "chr" prefix is optional and auto-matched to the BAM.
REGION_CHR = "9"
REGION_START = 89185060
REGION_END = 90828384
REGION_NAME = "Rasgrf1"

CONTROL_CONTIG_RE = re.compile(
    r"^(phage_lambda|plasmid_puc19c|Lambda|pUC19|J02459\.1)$"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "bismark_dir", help="Directory holding <sample>/<sample>.deduplicated.bam"
    )
    parser.add_argument(
        "--ref-fasta",
        default="",
        metavar="FASTA",
        help="Full reference to carve the mini genome from (e.g. "
        "resources/combined_genome.fa). Omit to skip mini-genome creation.",
    )
    parser.add_argument(
        "--ref-gtf",
        default="",
        metavar="GTF",
        help="Optional GTF (plain or .gz) to carve a matching mini annotation from. "
        "Its coordinates get the same shift/rename as the mini genome.",
    )
    parser.add_argument(
        "--ref-regulatory-gff3",
        default="",
        metavar="GFF3",
        help="Optional Ensembl regulatory-build GFF3 (plain or .gz) to carve a "
        "matching mini regulatory-features file from. Feeds generate_regions.py's "
        'promoter extraction (see resources.py\'s "test" genome branch).',
    )
    parser.add_argument(
        "--ref-cpg-islands-table",
        default="",
        metavar="TXT",
        help="Optional UCSC goldenPath cpgIslandExt.txt(.gz) table dump to carve a "
        "matching mini CpG islands track from.",
    )
    parser.add_argument(
        "--ref-repeat-mask-table",
        default="",
        metavar="TXT",
        help="Optional UCSC goldenPath rmsk.txt(.gz) table dump to carve a matching "
        "mini RepeatMasker track from.",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=100000,
        metavar="BP",
        help="bp added on each side of the target region (default: %(default)s)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads for samtools/pysam decompression (default: %(default)s)",
    )
    parser.add_argument(
        "--no-controls",
        action="store_true",
        help="Skip spike-in control reads entirely (default: include, capped at "
        "--control-max). Controls are kept so the methylation_conversion_rate rule "
        "has input.",
    )
    parser.add_argument(
        "--control-max",
        type=int,
        default=2000,
        metavar="N",
        help="Cap kept control reads per sample (default: %(default)s)",
    )
    parser.add_argument(
        "--subsample",
        default="",
        metavar="SEED.FRAC",
        help="Extra downsampling of the *region* reads, passed to `samtools view "
        "-s` (e.g. 42.25 keeps 25%%). Mates stay paired. Default: keep all.",
    )
    parser.add_argument(
        "--keep-coords",
        action="store_true",
        help="N-pad the 5' side so the mini genome keeps real mm39 coordinates "
        "(bigger file, but TxDb annotation / real BED files line up). Default: "
        "compact contig, region shifted to start near position 1.",
    )
    args = parser.parse_args()

    args.mini_contig = REGION_CHR.removeprefix("chr")
    args.win_start = max(REGION_START - args.flank, 1)
    args.win_end = REGION_END + args.flank
    # shift applied to every coordinate written into the test data
    args.offset = 0 if args.keep_coords else args.win_start - 1
    return args


def run(*cmd):
    subprocess.run(cmd, check=True)


def run_piped(cmds):
    """Runs `cmds` (a list of argv lists) as a shell-less unix pipeline."""
    procs = []
    stdin = None
    for cmd in cmds:
        p = subprocess.Popen(cmd, stdin=stdin, stdout=subprocess.PIPE)
        if stdin is not None:
            stdin.close()
        stdin = p.stdout
        procs.append(p)
    for p in procs:
        p.wait()
    for cmd, p in zip(cmds, procs):
        if p.returncode != 0:
            raise subprocess.CalledProcessError(p.returncode, cmd)


def resolve_contig(candidates, available):
    for candidate in candidates:
        if candidate in available:
            return candidate
    return None


def contig_candidates(chrom):
    bare = chrom.removeprefix("chr")
    return (bare, f"chr{bare}", chrom)


# ============================ per-sample FASTQ ================================


def scan_bam(args, bam_path):
    """
    Single pass over every record: buckets QNAMEs into "overlaps the target
    window" vs. "on a spike-in control contig" (mirroring what a
    `samtools view | awk` scan would do, but reading binary records
    directly instead of round-tripping through SAM text), and reports
    whether the file is paired-end from the first primary record seen.
    """
    with pysam.AlignmentFile(bam_path, "rb", threads=args.threads) as bam:
        contig = resolve_contig(contig_candidates(REGION_CHR), set(bam.references))
        if contig is None:
            sys.exit(f"ERROR: no contig for '{REGION_CHR}' in {bam_path}")
        # Matching the control regex against every reference name up front
        # (there are only ever a handful of contigs) turns the per-record
        # check into a cheap set lookup instead of a regex call per read.
        control_refs = {r for r in bam.references if CONTROL_CONTIG_RE.match(r)}

        region_names, control_names = set(), set()
        paired = None
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            if paired is None:
                paired = read.is_paired

            if read.is_unmapped:
                continue
            ref_name = read.reference_name
            if ref_name == contig:
                pos = read.reference_start + 1  # 1-based, matches SAM POS
                length = read.query_length or 0
                if pos <= args.win_end and (pos + length) >= args.win_start:
                    region_names.add(read.query_name)
                    continue
            if ref_name in control_refs:
                control_names.add(read.query_name)

    return contig, region_names, control_names, bool(paired)


def extract_fastq(args, bam_path, names, paired, out_prefix):
    with tempfile.NamedTemporaryFile("wt", suffix=".txt", delete=False) as fh:
        fh.write("\n".join(sorted(names)) + "\n")
        names_file = fh.name

    threads = str(args.threads)
    view_cmd = [
        "samtools",
        "view",
        "-@",
        threads,
        "-b",
        "-F",
        "0x900",
        "-N",
        names_file,
        *(["-s", args.subsample] if args.subsample else []),
        str(bam_path),
    ]
    try:
        if paired:
            run_piped(
                [
                    view_cmd,
                    ["samtools", "collate", "-@", threads, "-u", "-O", "-"],
                    [
                        "samtools",
                        "fastq",
                        "-@",
                        threads,
                        "-n",
                        "-1",
                        f"{out_prefix}_R1_001.fastq.gz",
                        "-2",
                        f"{out_prefix}_R2_001.fastq.gz",
                        "-0",
                        "/dev/null",
                        "-s",
                        "/dev/null",
                        "-",
                    ],
                ]
            )
        else:
            run_piped(
                [
                    view_cmd,
                    [
                        "samtools",
                        "fastq",
                        "-@",
                        threads,
                        "-n",
                        "-0",
                        f"{out_prefix}.fastq.gz",
                        "-",
                    ],
                ]
            )
    finally:
        os.remove(names_file)


def make_reads(args):
    bams = sorted(glob.glob(str(Path(args.bismark_dir) / "*" / "*.deduplicated.bam")))
    if not bams:
        sys.exit(f"ERROR: no *.deduplicated.bam under {args.bismark_dir}")

    READS_DIR.mkdir(parents=True, exist_ok=True)
    for old in READS_DIR.glob("*.fastq.gz"):
        old.unlink()

    print(
        f"Window +flank: {REGION_CHR}:{args.win_start}-{args.win_end}  "
        f"({args.flank} bp each side)"
    )
    print(f"Output       : {READS_DIR}\n")

    for bam_path in bams:
        sample = Path(bam_path).parent.name
        print(f">> {sample}")

        contig, region_names, control_names, paired = scan_bam(args, bam_path)
        kept_controls = (
            [] if args.no_controls else sorted(control_names)[: args.control_max]
        )
        names = region_names | set(kept_controls)
        print(
            f"   region read names: {len(region_names)}   "
            f"control read names: {len(kept_controls)}"
        )
        if not region_names:
            print(f"   WARNING: no reads found in region for {sample}")

        extract_fastq(args, bam_path, names, paired, str(READS_DIR / sample))


# ============================ mini genome ====================================


def make_mini_genome(args):
    ref_fasta = args.ref_fasta
    if not ref_fasta:
        print(">> --ref-fasta not set -> skipped mini genome.")
        print("   Re-run with e.g.  --ref-fasta resources/combined_genome.fa")
        return
    if not Path(ref_fasta).is_file():
        print(f">> --ref-fasta '{ref_fasta}' not found -> skipped mini genome.")
        return

    print(f">> mini genome from {ref_fasta}")
    if not Path(f"{ref_fasta}.fai").exists():
        pysam.faidx(ref_fasta)

    with pysam.FastaFile(ref_fasta) as fasta:
        rcontig = resolve_contig(contig_candidates(REGION_CHR), set(fasta.references))
        if rcontig is None:
            sys.exit(f"   ERROR: '{REGION_CHR}' not found in {ref_fasta}.fai")
        seq = fasta.fetch(rcontig, args.win_start - 1, args.win_end)

    if args.keep_coords:
        seq = "N" * (args.win_start - 1) + seq

    out_fasta = RES_DIR / "genome.fa"
    with open(out_fasta, "wt") as out:
        out.write(f">{args.mini_contig}\n")
        out.write("\n".join(textwrap.wrap(seq, 60)))
        out.write("\n")
    pysam.faidx(str(out_fasta))

    # resources.py's "test" genome branch fetches genome.fa.gz (matching
    # every other genome's fasta_url, which always points at a .gz) -- keep
    # this in sync with the plain genome.fa written above, the same way
    # carve_annotation() below does for the annotation tracks.
    with open(out_fasta, "rb") as f_in, gzip.open(f"{out_fasta}.gz", "wb") as f_out:
        f_out.writelines(f_in)

    chrom_sizes = RES_DIR / "chrom.sizes"
    with open(f"{out_fasta}.fai") as fai, open(chrom_sizes, "wt") as out:
        for line in fai:
            name, length = line.split("\t")[:2]
            out.write(f"{name}\t{length}\n")

    target_bed = RES_DIR / "target_region.bed"
    with open(target_bed, "wt") as out:
        out.write(
            f"{args.mini_contig}\t{REGION_START - args.offset - 1}\t"
            f"{REGION_END - args.offset}\t{REGION_NAME}\n"
        )

    print(
        f"   {out_fasta}(.gz)  (contig '{args.mini_contig}', coordinate offset "
        f"applied to test data = {args.offset})"
    )
    print(f"   {target_bed}  (target locus in mini-genome coordinates)")


# ================= mini annotation tracks (GTF/GFF3/UCSC tables) =============


def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def carve_annotation(args, label, ref_path, out_path, chrom_col, start_col, end_col):
    """
    Keeps rows overlapping the target window, clamps overhangs, shifts +
    renames the chrom column to match the mini genome. Column positions
    (0-based) differ per format:
        GTF/GFF3       : chrom_col=0, start_col=3, end_col=4
        cpgIslandExt   : chrom_col=1, start_col=2, end_col=3
        rmsk           : chrom_col=5, start_col=6, end_col=7
    All other columns are passed through unchanged.
    """
    if not ref_path:
        return
    if not Path(ref_path).is_file():
        print(f">> {label}: '{ref_path}' not found -> skipped.")
        return
    print(f">> {label} from {ref_path}")

    candidates = set(contig_candidates(REGION_CHR))
    lo, hi = args.win_start - args.offset, args.win_end - args.offset
    contig = None
    n = 0
    with open_maybe_gzip(ref_path) as fh, open(out_path, "wt") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[chrom_col]
            if contig is None:
                if chrom not in candidates:
                    continue
                contig = chrom
            elif chrom != contig:
                continue

            start, end = int(fields[start_col]), int(fields[end_col])
            if start > args.win_end or end < args.win_start:
                continue

            fields[chrom_col] = args.mini_contig
            fields[start_col] = str(max(start - args.offset, lo))
            fields[end_col] = str(min(end - args.offset, hi))
            out.write("\t".join(fields) + "\n")
            n += 1

    if contig is None:
        sys.exit(f"   ERROR: no contig for '{REGION_CHR}' in {ref_path}")

    with open(out_path, "rb") as f_in, gzip.open(f"{out_path}.gz", "wb") as f_out:
        f_out.writelines(f_in)
    print(
        f"   {out_path}(.gz)  ({n} features overlapping the window, overhangs clamped)"
    )


def make_annotations(args):
    carve_annotation(args, "mini GTF", args.ref_gtf, RES_DIR / "genes.gtf", 0, 3, 4)
    carve_annotation(
        args,
        "mini regulatory GFF3",
        args.ref_regulatory_gff3,
        RES_DIR / "regulatory_features.gff3",
        0,
        3,
        4,
    )
    carve_annotation(
        args,
        "mini CpG islands table",
        args.ref_cpg_islands_table,
        RES_DIR / "cpgIslandExt.txt",
        1,
        2,
        3,
    )
    carve_annotation(
        args,
        "mini RepeatMasker table",
        args.ref_repeat_mask_table,
        RES_DIR / "rmsk.txt",
        5,
        6,
        7,
    )


def main():
    args = parse_args()

    RES_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Region      : {REGION_CHR}:{REGION_START}-{REGION_END} ({REGION_NAME})")
    make_reads(args)

    print()
    make_mini_genome(args)
    print()
    make_annotations(args)

    print("\nDone.")
    outputs = sorted(str(p) for p in READS_DIR.glob("*"))
    if outputs:
        run("du", "-sh", *outputs)


if __name__ == "__main__":
    main()
