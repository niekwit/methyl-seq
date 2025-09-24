#!/usr/bin/env python3

"""
This script identifies CpG dinucleotide sites ("CG") in a given genome FASTA file and outputs their locations in BED format (gzipped).
It processes each chromosome in parallel using multiprocessing for efficiency.

Usage:
    python find_cpgs.py <input_genome.fasta> <output.bed> <log_file>

Arguments:
    input_genome.fasta : Path to the input genome FASTA file.
    output.bed         : Path to the output gzipped BED file containing CpG locations.

Output:
    A gzipped BED file listing all CpG sites in the genome, with each line containing:
        <chromosome_name> <start_position> <end_position>
    Positions are 0-based and end-exclusive, following BED format conventions.

Dependencies:
    - pyfaidx
    - tqdm
    - gzip
    - multiprocessing
    - sys

Functions:
    process_chromosome(args):
        Finds all CpG sites in a given chromosome and returns their BED-formatted locations.

Notes:
    - The script uses multiprocessing to process chromosomes in parallel.
    - Progress is displayed using tqdm.
    - The output file is compressed using gzip.
"""

from pyfaidx import Fasta
import sys
import multiprocessing
import logging


def process_chromosome(args):
    """
    Worker function to find CpG sites in a single chromosome.
    Returns a list of strings, one for each CpG location.
    """
    chrom_name, fasta_file = args
    cpg_locations = []

    # Open the FASTA file. Pyfaidx is thread-safe for reading.
    try:
        genome = Fasta(fasta_file)
        chrom_seq = genome[chrom_name]

        # Iterate through the sequence to find "CG" dinucleotides.
        logging.info(f"Processing chromosome: {chrom_name}")
        for i in range(len(chrom_seq) - 1):
            if str(chrom_seq[i : i + 2]).upper() == "CG":
                # BED format is 0-based, so the start position is 'i'
                # The end position is 'i+2' (exclusive)
                cpg_locations.append(f"{chrom_name}\t{i}\t{i+2}\n")
    except Exception as e:
        logging.error(f"Error processing chromosome {chrom_name}: {e}")
        return []

    return cpg_locations


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python find_cpgs.py <input_genome.fasta.gz> <output.bed> <log_file>"
        )
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_bed = sys.argv[2]
    log = sys.argv[3]

    # Set up logging
    logging.basicConfig(
        format="%(levelname)s:%(asctime)s: %(message)s",
        level=logging.DEBUG,
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.FileHandler(log)],
    )

    logging.info("Starting chromosome processing...")
    logging.info(f"Input FASTA file: {input_fasta}")
    logging.info(f"Output BED file: {output_bed}")

    try:
        # Get all chromosome names. This is done once to avoid parallel issues.
        all_chroms = list(Fasta(input_fasta).keys())
    except FileNotFoundError:
        logging.error(f"The file '{input_fasta}' was not found.")
        sys.exit(1)

    # Set the number of processes to use.
    num_processes = 10

    # Create a pool of worker processes.
    with multiprocessing.Pool(processes=num_processes) as pool:
        args_list = [(chrom, input_fasta) for chrom in all_chroms]
        all_results = pool.map(process_chromosome, args_list)

    # Write all results to the output file.
    logging.info("Finished processing all chromosomes")
    logging.info(f"Writing results to {output_bed}")
    with open(output_bed, "wt") as bed_file:
        for result_list in all_results:
            bed_file.writelines(result_list)

    logging.info(f"Successfully wrote all CpG locations to file!")
