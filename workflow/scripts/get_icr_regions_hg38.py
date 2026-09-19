#!/usr/bin/env python3

"""
Fetch the 25 canonical human imprinting control regions (ICRs) in hg38
coordinates, fully referenced, and write them as a BED file plus a
companion references table.

Data sources (queried live, not hardcoded):

  1. Coordinates + region IDs: the "known_ICR_25" track of the Human ICRs
     database (https://humanicr.org/), a JBrowse1-backed resource built
     for:
         Sanchez-Delgado M, Court F, Vidal-Bralo L, et al. "Genomic map
         of candidate human imprint control regions: the imprintome."
         Epigenetics. 2022;17(13):2181-2197. PMID: 35786392.
     This paper explicitly built known_ICR_25 from the pre-existing,
     literature-curated set of 25 characterized human ICRs (see #2) --
     it is not this paper's own novel candidate list (that is the
     separate, much larger "putative_ICRs_35_65_final_N1488" track).
     Reference genome build: GRCh38 (confirmed here by chr1 length ==
     248,956,422 bp, the canonical GRCh38 chr1 length).

  2. The original enumeration of these 25 regions as "known"/characterized
     human ICRs:
         Skaar DA, Li Y, Bernal AJ, Hoyo C, Murphy SK, Jirtle RL. "The
         human imprintome: regulatory mechanisms, methods of
         ascertainment, and roles in disease susceptibility." ILAR J.
         2012;53(3-4):341-358. PMID: 23744971.

  3. Gene-symbol annotation for each region (the humanicr.org track only
     carries opaque IDs like "Known_ICR_9", not gene names): the UCSC
     Genome Browser REST API (https://api.genome.ucsc.edu), querying the
     hg38 ncbiRefSeq track for the gene(s) overlapping each ICR interval.

Output:
  --out        BED4 (chrom, start, end, name) of all 25 regions, name =
               the overlapping gene symbol/locus (falling back to the
               source database's own "Known_ICR_N" ID if no NCBI RefSeq
               gene overlaps that exact interval).
  --refs-out   TSV with one row per region carrying its coordinates, ID,
               resolved gene symbol(s), and the full reference list above
               -- so every row is traceable back to its source without
               needing this script's docstring.

This is a standalone, manually-run utility (not part of the Snakemake
DAG), the same way .test/make_test_data.py is -- it hits two live,
third-party web services, which is inappropriate to do automatically as
part of a pipeline run.
"""

import argparse
import logging
import time
import urllib.error
import urllib.request
import json

HUMANICR_BASE = "https://jb2.humanicr.org/data/json/humanicr"
ICR_TRACK = "known_ICR_25"
UCSC_API = "https://api.genome.ucsc.edu/getData/track"
UCSC_GENOME = "hg38"
UCSC_GENE_TRACK = "ncbiRefSeq"

# GRCh38 primary assembly chromosomes only -- known_ICR_25 is a small,
# hand-curated set, so there is no reason to expect it to place anything
# on an alt/patch/scaffold contig, and querying only these keeps this
# script's runtime and request count small.
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

REFERENCES = {
    "coordinates_and_list": (
        "Sanchez-Delgado M, Court F, Vidal-Bralo L, et al. Genomic map of "
        "candidate human imprint control regions: the imprintome. "
        "Epigenetics. 2022;17(13):2181-2197. PMID: 35786392. "
        "Data: https://humanicr.org/ (known_ICR_25 track)."
    ),
    "original_25_icr_characterization": (
        "Skaar DA, Li Y, Bernal AJ, Hoyo C, Murphy SK, Jirtle RL. The human "
        "imprintome: regulatory mechanisms, methods of ascertainment, and "
        "roles in disease susceptibility. ILAR J. 2012;53(3-4):341-358. "
        "PMID: 23744971."
    ),
    "gene_symbol_annotation": (
        "UCSC Genome Browser REST API, hg38 ncbiRefSeq track. "
        "https://api.genome.ucsc.edu"
    ),
}


def fetch_json(url, retries=3, pause=1.0):
    last_err = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                return json.loads(resp.read())
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            last_err = e
            logging.debug(f"Fetch failed ({attempt + 1}/{retries}) for {url}: {e}")
            time.sleep(pause)
    if last_err is not None:
        raise last_err


def assert_hg38(chrom1_len):
    if chrom1_len != 248956422:
        raise SystemExit(
            f"ERROR: humanicr.org's chr1 length is {chrom1_len} bp, not "
            "GRCh38's 248,956,422 bp -- the database's reference build may "
            "have changed. Refusing to emit possibly-non-hg38 coordinates; "
            "check https://humanicr.org/ manually before re-running."
        )


def fetch_known_icrs():
    """Query known_ICR_25 across all standard chromosomes.

    Returns a list of dicts: {chrom, start, end, source_id} in hg38
    0-based half-open BED coordinates.
    """
    refseqs = fetch_json(f"{HUMANICR_BASE}/seq/refSeqs.json")
    chrom1 = next((r for r in refseqs if r["name"] == "chr1"), None)
    if chrom1 is None:
        raise SystemExit("ERROR: chr1 not found in humanicr.org's refSeqs.json")
    assert_hg38(chrom1["length"])

    regions = []
    for chrom in CHROMS:
        url = f"{HUMANICR_BASE}/tracks/{ICR_TRACK}/{chrom}/trackData.json"
        try:
            data = fetch_json(url)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                # No known_ICR_25 features on this chromosome.
                continue
            raise
        if not isinstance(data, dict) or data.get("featureCount", 0) == 0:
            continue
        # NCList row layout per trackList.json's declared attributes:
        # [level, Start, End, Strand, Name, Seq_id]
        for row in data["intervals"]["nclist"]:
            start, end, name = row[1], row[2], row[4]
            regions.append(
                {"chrom": chrom, "start": start, "end": end, "source_id": name}
            )
        logging.info(f"{chrom}: {data['featureCount']} known ICR(s)")

    regions.sort(key=lambda r: (CHROMS.index(r["chrom"]), r["start"]))
    return regions


def ensembl_chrom(chrom):
    """UCSC/humanicr.org use "chr1"-style names; this pipeline's other hg38
    resources (fasta, GTF, repeat_mask.bed, cpg_islands.bed -- all fetched
    from Ensembl or stripped of "chr" from a UCSC source, see
    prepare_repeat_mask/prepare_cpg_islands in resources.smk) use bare
    Ensembl-style names ("1", "X", "MT"). Match that convention so this
    file is directly usable alongside them."""
    name = chrom[3:] if chrom.startswith("chr") else chrom
    return "MT" if name == "M" else name


def resolve_gene_symbol(chrom, start, end):
    """Return a gene symbol overlapping [start, end) on chrom via the UCSC
    ncbiRefSeq track, or None if nothing overlaps."""
    url = (
        f"{UCSC_API}?genome={UCSC_GENOME};track={UCSC_GENE_TRACK};"
        f"chrom={chrom};start={start};end={end}"
    )
    try:
        data = fetch_json(url)
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        logging.warning(f"UCSC lookup failed for {chrom}:{start}-{end}: {e}")
        return None
    features = data.get(UCSC_GENE_TRACK, [])
    symbols = []
    for f in features:
        sym = f.get("name2")
        if sym and sym not in symbols:
            symbols.append(sym)
    if not symbols:
        return None
    return "/".join(symbols)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--out", required=True, help="Output BED4 file (chrom, start, end, name)"
    )
    parser.add_argument(
        "--refs-out",
        required=True,
        help="Output TSV with per-region coordinates, source ID, resolved "
        "gene symbol, and full references",
    )
    parser.add_argument(
        "--no-gene-lookup",
        action="store_true",
        help="Skip the UCSC gene-symbol lookup and name regions by their "
        "humanicr.org source_id only (faster, no second web service)",
    )
    parser.add_argument("--log", help="Log file (default: stderr)")
    return parser.parse_args()


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

    logging.info("Fetching known_ICR_25 from humanicr.org ...")
    regions = fetch_known_icrs()
    if len(regions) != 25:
        logging.warning(
            f"Expected 25 known ICRs, got {len(regions)} -- humanicr.org's "
            "content may have changed since this script was written."
        )
    logging.info(f"Fetched {len(regions)} region(s)")

    if not args.no_gene_lookup:
        logging.info("Resolving gene symbols via UCSC REST API ...")
        for r in regions:
            r["gene"] = resolve_gene_symbol(r["chrom"], r["start"], r["end"])
            time.sleep(0.2)  # be polite to the public API
    else:
        for r in regions:
            r["gene"] = None

    with open(args.out, "w") as bed:
        for r in regions:
            name = r["gene"] or r["source_id"]
            bed.write(
                f"{ensembl_chrom(r['chrom'])}\t{r['start']}\t{r['end']}\t{name}\n"
            )
    logging.info(f"Wrote {len(regions)} region(s) to {args.out}")

    with open(args.refs_out, "w") as tsv:
        tsv.write(
            "chrom\tstart\tend\tsource_id\tgene_symbol\t"
            "reference_coordinates_and_list\treference_original_25_icr_list\t"
            "reference_gene_annotation\n"
        )
        for r in regions:
            tsv.write(
                f"{ensembl_chrom(r['chrom'])}\t{r['start']}\t{r['end']}\t{r['source_id']}\t"
                f"{r['gene'] or ''}\t"
                f"{REFERENCES['coordinates_and_list']}\t"
                f"{REFERENCES['original_25_icr_characterization']}\t"
                f"{REFERENCES['gene_symbol_annotation']}\n"
            )
    logging.info(f"Wrote references table to {args.refs_out}")


if __name__ == "__main__":
    main()
