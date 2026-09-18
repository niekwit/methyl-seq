#!/usr/bin/env python3

"""
Locate ORF1/ORF2 and define 5'/3' UTR boundaries for near-full-length LINE1
(L1) elements.

Each candidate L1 element (from --bed) is treated as a candidate full-length
copy. For each one, the ORF1 and ORF2 consensus nucleotide sequences (read
from a GenBank record) are mapped onto that element's own genomic sequence
with minimap2 (via mappy) to find where ORF1/ORF2 actually sit in that
specific copy. The 5' UTR is defined as the sequence upstream of ORF1 and the
3' UTR as the sequence downstream of ORF2, both in the element's own 5'->3'
sense (i.e. strand-aware).

A fifth feature, "remainder", is also written: the rest of the element that
is not the 5' UTR (i.e. ORF1 + linker + ORF2 + 3' UTR as one span).

Elements where ORF1 and/or ORF2 cannot be confidently located (low identity,
low coverage, ambiguous mapping, or an ORF1/ORF2 order inconsistent with the
annotated strand) are skipped and reported in the summary/log. Since ORF
localisation success is data-dependent (sequence divergence, element
completeness), zero annotated elements is not treated as an error here --
unlike a class/family/subfamily filter matching nothing, it doesn't
necessarily indicate a config mistake.

An ORF hit must clear two independent thresholds:
  --min-identity = matched bases / aligned block length: how similar the copy
    is to the consensus *within the region that aligned*.
  --min-coverage = aligned block length / full consensus ORF length: how much
    of the consensus ORF aligned at all.
They catch different bad hits: a truncated ORF can align at high identity but
low coverage, and a full-length but heavily mutated ORF can have high coverage
but low identity. Low-coverage hits are especially dangerous here because they
place the ORF end - and hence the UTR boundary - in the wrong spot.

By default an element needs BOTH ORFs to be annotated. With --single-orf,
elements where exactly one of ORF1/ORF2 is confidently located are still
annotated - from that single ORF alone, with the 5'/3' UTR defined as the
sequence up-/downstream of it - and written to --out alongside the
fully-annotated elements instead of being skipped. Such an element just
contributes fewer rows (none for the unlocated ORF). Because the other ORF
is unlocated, the "UTR" on its side also contains the missing ORF and the
inter-ORF sequence.

Suggestions for LINE1 GenBank files:
Mouse: https://www.ncbi.nlm.nih.gov/nuccore/M13002
Human: https://www.ncbi.nlm.nih.gov/nuccore/AF148856
Note: /product qualifiers must be exactly "ORF1" and "ORF2" for the script to find them.
"""

import argparse
import logging
import re
import tempfile
from pathlib import Path

import mappy
import pyfaidx
from Bio import SeqIO

SCAFFOLD_RE = re.compile(r"^(?:[0-9XY]+_|Un_)([A-Za-z]+\d+)v(\d+)(?:_random)?$")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", required=True, help="BED6+ file of LINE1 elements")
    parser.add_argument(
        "--genbank",
        required=True,
        help="GenBank record with ORF1 and ORF2 CDS features (consensus L1)",
    )
    parser.add_argument("--fasta", required=True, help="Genome FASTA (faidx-indexed)")
    parser.add_argument("--out", required=True, help="Output BED6 file")
    parser.add_argument(
        "--min-length",
        type=int,
        default=6000,
        help="Minimum element length (bp) to consider (default: 6000)",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=20,
        help="Minimum mappy mapping quality for an ORF hit (default: 20)",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.6,
        help="Minimum identity (matched/aligned bases) for an ORF hit (default: 0.6)",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.5,
        help="Minimum fraction of the consensus ORF length covered by the "
        "hit (default: 0.5)",
    )
    parser.add_argument(
        "--skipped-out",
        help="Optional TSV path to write skipped elements and the reason",
    )
    parser.add_argument(
        "--single-orf",
        action="store_true",
        help="Also annotate elements where only ONE of ORF1/ORF2 is "
        "confidently located, from that ORF alone (5'/3' UTR = sequence "
        "up-/downstream of it), writing them to --out alongside the "
        "fully-annotated elements instead of skipping them. Such an element "
        "contributes fewer rows (none for the unlocated ORF), and the UTR on "
        "the side of the missing ORF also contains that ORF and the inter-ORF "
        "sequence.",
    )
    parser.add_argument("--log", help="Log file (default: stderr)")
    return parser.parse_args()


def load_orf_consensus(genbank_path):
    """Return {"ORF1": seq, "ORF2": seq} from the consensus L1 GenBank record."""
    record = SeqIO.read(genbank_path, "genbank")
    orfs = {}
    # Walk every annotated feature and keep the two CDS entries whose /product
    # qualifier is exactly "ORF1" or "ORF2", storing their spliced nucleotide
    # sequence.
    for feature in record.features:
        if feature.type != "CDS":
            continue
        product = feature.qualifiers.get("product", [None])[0]
        if product not in ("ORF1", "ORF2"):
            continue
        orfs[product] = str(feature.location.extract(record.seq))
    # Both ORFs are mandatory; a record missing either one cannot drive the
    # downstream mapping, so fail loudly here.
    missing = {"ORF1", "ORF2"} - orfs.keys()
    if missing:
        raise ValueError(
            f"{genbank_path} is missing CDS feature(s) with product "
            f"{sorted(missing)}"
        )
    return orfs


def build_aligner(orfs, tmp_dir):
    """Index the ORF consensus sequences as minimap2 references (mappy).

    The consensus ORFs become the *reference* and each element's own genomic
    sequence is later mapped as a *query* against them, so a hit tells us where
    ORF1/ORF2 sit within that element.
    """
    # mappy indexes from a file on disk, so dump the consensus sequences to a
    # throwaway FASTA in the caller's temp dir.
    fasta_path = Path(tmp_dir) / "orf_consensus.fa"
    with open(fasta_path, "w") as fh:
        for name, seq in orfs.items():
            fh.write(f">{name}\n{seq}\n")
    # "map-ont" is a permissive long-read preset; it tolerates the substantial
    # divergence between the consensus and individual (often old, mutated) L1
    # copies better than the short-read presets.
    aligner = mappy.Aligner(str(fasta_path), preset="map-ont")
    if not aligner:
        raise RuntimeError(
            "Failed to build minimap2 index from ORF consensus sequences"
        )
    return aligner


def normalize_chrom(bed_chrom, fasta_chroms):
    """Map a RepeatMasker/UCSC-style scaffold name to its Ensembl FASTA name.

    Returns the matching FASTA contig name, or None if the BED contig has no
    counterpart in the genome FASTA (caller then skips that element).
    """
    # Common case: the BED already uses the same names as the FASTA.
    if bed_chrom in fasta_chroms:
        return bed_chrom
    # Otherwise try to rewrite a UCSC-style unplaced-scaffold name such as
    # "1_GL456210v1_random" or "Un_GL456210v1" -> Ensembl "GL456210.1".
    match = SCAFFOLD_RE.match(bed_chrom)
    if match:
        candidate = f"{match.group(1)}.{match.group(2)}"
        if candidate in fasta_chroms:
            return candidate
    return None


def best_orf_hits(aligner, seq, orf_lengths, min_mapq, min_identity, min_coverage):
    """Map one element's sequence to the ORF consensus and keep the best hit
    per ORF that clears all quality thresholds.

    Returns {"ORF1": hit, "ORF2": hit} with 0, 1 or 2 entries.
    """
    best = {}
    # aligner.map yields every alignment minimap2 finds between this element and
    # the two consensus ORFs; filter them one ORF at a time.
    for hit in aligner.map(seq):
        # Reference name must be one of our ORFs (guards against stray contigs).
        if hit.ctg not in orf_lengths:
            continue
        # Drop low-confidence / multi-mapping placements.
        if hit.mapq < min_mapq:
            continue
        # Identity = matching bases / aligned block length.
        if hit.mlen / hit.blen < min_identity:
            continue
        # Coverage = aligned block length / full consensus ORF length; rejects
        # short partial hits that would misplace the UTR boundary.
        if hit.blen / orf_lengths[hit.ctg] < min_coverage:
            continue
        # Among the survivors for this ORF, keep the one with the most matched
        # bases.
        current = best.get(hit.ctg)
        if current is None or hit.mlen > current.mlen:
            best[hit.ctg] = hit
    return best


def _orf_span(start, orf_hit):
    """Lift a hit's query offsets (into the element sequence) to genomic coords."""
    return (start + orf_hit.q_st, start + orf_hit.q_en)


def _features_from_anchor(start, end, strand, upstream_end, downstream_start, spans):
    """Assemble the feature dict shared by the two- and single-ORF cases.

    `spans` is {ORF name: (gstart, gend)}; `upstream_end`/`downstream_start` are
    the genomic coordinates that bound the 5' UTR (upstream, toward the element's
    5' end) and the 3' UTR (downstream, toward its 3' end).
    """
    if strand == "+":
        # Forward element: 5' end is `start`, 3' end is `end`.
        utr5 = (start, upstream_end)
        features = {"5UTR": utr5, **spans, "3UTR": (downstream_start, end)}
        features["remainder"] = (utr5[1], end)
    else:
        # Reverse element: 5' end is `end`, 3' end is `start`.
        utr5 = (downstream_start, end)
        features = {"5UTR": utr5, **spans, "3UTR": (start, upstream_end)}
        features["remainder"] = (start, utr5[0])
    return features


def locate_features(start, end, strand, hits):
    """Turn ORF1/ORF2 hits into genomic coordinates for the element's features.

    `start`/`end` are the element's genomic coordinates; `hits` holds the ORF
    positions as offsets into the element sequence (q_st/q_en). Returns
    (features, used_orfs, None) on success - where `features` maps a feature
    name to (genomic_start, genomic_end) and `used_orfs` is {"ORF1", "ORF2"}
    for a full annotation or a single-element set when only one ORF was found -
    or (None, None, reason) if the element fails a sanity check.
    """
    # mappy strand: +1 = forward, -1 = reverse; must match the annotated strand.
    expected_strand = 1 if strand == "+" else -1
    present = [name for name in ("ORF1", "ORF2") if name in hits]
    if not present:
        return None, None, "ORF1 or ORF2 not confidently located"

    # Every located ORF must map in the element's annotated orientation.
    for name in present:
        if hits[name].strand != expected_strand:
            return None, None, "ORF hit strand disagrees with annotated element strand"

    if len(present) == 2:
        orf1 = _orf_span(start, hits["ORF1"])
        orf2 = _orf_span(start, hits["ORF2"])
        # Along increasing genomic coordinates the coding pair reads ORF1->ORF2
        # on a + element and ORF2->ORF1 on a - element; anything else means the
        # two hits are inconsistent and the element is dropped.
        if strand == "+":
            if orf1[1] > orf2[0]:
                return None, None, "ORF1/ORF2 out of expected 5'->3' order"
            upstream_end, downstream_start = orf1[0], orf2[1]
        else:
            if orf2[1] > orf1[0]:
                return None, None, "ORF1/ORF2 out of expected 5'->3' order"
            upstream_end, downstream_start = orf2[0], orf1[1]
        spans = {"ORF1": orf1, "ORF2": orf2}
        features = _features_from_anchor(
            start, end, strand, upstream_end, downstream_start, spans
        )
        return features, frozenset({"ORF1", "ORF2"}), None

    # Single-ORF fallback: exactly one ORF located. The UTRs are simply the
    # sequence up- and downstream of that ORF, so the "UTR" on the side of the
    # missing ORF also swallows that ORF and the inter-ORF sequence.
    name = present[0]
    orf = _orf_span(start, hits[name])
    features = _features_from_anchor(start, end, strand, orf[0], orf[1], {name: orf})
    return features, frozenset({name}), None


def main():
    args = parse_args()

    logging.basicConfig(
        format="%(levelname)s:%(asctime)s: %(message)s",
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=(
            [logging.FileHandler(args.log)] if args.log else [logging.StreamHandler()]
        ),
    )
    logging.info(f"Arguments: {vars(args)}")

    # Load the consensus ORF sequences and record their lengths (needed for the
    # coverage filter later).
    orfs = load_orf_consensus(args.genbank)
    orf_lengths = {name: len(seq) for name, seq in orfs.items()}
    lengths_str = {k: f"{v} bp" for k, v in orf_lengths.items()}
    logging.info(f"Loaded ORF consensus: {lengths_str}")

    # Open the genome once; keep the set of contig names for name normalisation.
    genome = pyfaidx.Fasta(args.fasta)
    fasta_chroms = set(genome.keys())

    n_total = n_candidates = n_annotated = n_single_orf = n_skipped = 0
    skipped_rows = []

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Build the minimap2 index of the two ORF consensus sequences once and
        # reuse it for every element.
        aligner = build_aligner(orfs, tmp_dir)

        with open(args.bed) as bed_in, open(args.out, "w") as bed_out:
            for line in bed_in:
                # Parse the first 6 BED columns; ignore any extra ones.
                fields = line.rstrip("\n").split("\t")
                chrom, start, end, name, _score, strand = fields[:6]
                start, end = int(start), int(end)
                n_total += 1

                # Step 1: only near-full-length copies are candidates.
                if end - start <= args.min_length:
                    continue
                n_candidates += 1

                # Step 2: resolve the BED contig name to a FASTA contig name.
                fasta_chrom = normalize_chrom(chrom, fasta_chroms)
                if fasta_chrom is None:
                    reason = "chromosome not in FASTA"
                    logging.info(f"{name} {chrom}:{start}-{end} skipped: {reason}")
                    n_skipped += 1
                    skipped_rows.append((chrom, start, end, name, reason))
                    continue

                # Step 3: pull this element's own genomic sequence and map the
                # consensus ORFs onto it.
                seq = str(genome[fasta_chrom][start:end])
                hits = best_orf_hits(
                    aligner,
                    seq,
                    orf_lengths,
                    args.min_mapq,
                    args.min_identity,
                    args.min_coverage,
                )
                # Step 4: convert ORF positions into feature intervals. `used`
                # is {"ORF1", "ORF2"} for a full annotation or a single-element
                # set when only one ORF was located.
                features, used, reason = locate_features(start, end, strand, hits)
                if features is None:
                    logging.info(f"{name} {chrom}:{start}-{end} skipped: {reason}")
                    n_skipped += 1
                    skipped_rows.append((chrom, start, end, name, reason))
                    continue

                # A single-ORF annotation is only kept when --single-orf was
                # given; otherwise keep the original all-or-nothing behaviour.
                single = len(used) < 2
                if single and not args.single_orf:
                    reason = f"only {next(iter(used))} confidently located"
                    logging.info(f"{name} {chrom}:{start}-{end} skipped: {reason}")
                    n_skipped += 1
                    skipped_rows.append((chrom, start, end, name, reason))
                    continue

                # Step 5: emit one BED row per feature, dropping any that came
                # out zero-length (e.g. an ORF flush against the element edge).
                # The feature name encodes the parent element and feature type.
                out_rows = [
                    (
                        chrom,
                        f_start,
                        f_end,
                        f"{name}_{chrom}_{start}_{end}_{feat}",
                        strand,
                    )
                    for feat, (f_start, f_end) in features.items()
                    if f_end > f_start
                ]
                # Write features in genomic-coordinate order; single-ORF
                # elements go to the same --out, just with fewer rows.
                out_rows.sort(key=lambda r: r[1])
                for c, s, e, n, st in out_rows:
                    bed_out.write(f"{c}\t{s}\t{e}\t{n}\t.\t{st}\n")
                if single:
                    n_single_orf += 1
                else:
                    n_annotated += 1

    if args.skipped_out and skipped_rows:
        with open(args.skipped_out, "w") as fh:
            fh.write("chrom\tstart\tend\tname\treason\n")
            for row in skipped_rows:
                fh.write("\t".join(str(x) for x in row) + "\n")

    logging.info(
        f"{n_total} elements total, {n_candidates} longer than {args.min_length} bp, "
        f"{n_annotated + n_single_orf} annotated ({n_single_orf} of them from a "
        f"single ORF), {n_skipped} skipped"
    )


if __name__ == "__main__":
    main()
