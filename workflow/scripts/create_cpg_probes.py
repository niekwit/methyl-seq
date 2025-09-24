#!/usr/bin/env python3


"""
This script generates non-overlapping genomic probes from a BED file containing CpG coordinates.
Each probe consists of a window of consecutive CpGs on the same chromosome, with user-defined window size.
Probes are created in parallel for each chromosome and written to a BED output file.

Usage:
    python create_cpg_probes.py <input.bed.gz> <output.bed.gz> <window_size> <log_file>

Arguments:
    input.bed.gz   Path to the BED file containing CpG coordinates (chrom, start, end).
    output.bed.gz  Path to the BED file to write the generated probes.
    window_size    Number of CpGs per probe (integer).
    log_file       Path to the log file for recording progress and errors.

Main Functions:
    - create_probes: Worker function to generate probes for a single chromosome.
    - get_chromosomes: Extracts unique chromosome names from the input BED file.

Features:
    - Parallel processing of chromosomes using multiprocessing.
    - Progress bar for probe creation using tqdm.
    - Sorted output by chromosome and start coordinate.
    - Robust error handling for file operations and data processing.

Output:
    A BED file where each line represents a probe with format:
    chrom    start    end    probe_name

Example probe line:
    1	3741327	3741745	probe_1_3741327_3741745
"""

import sys
import logging
import multiprocessing
from operator import itemgetter
import pandas as pd


def create_probes(args):
    """
    Worker function to create non-overlapping probes for a single chromosome.

    Args:
        args (tuple): A tuple containing the chromosome name, input file path,
                      and window size.

    Returns:
        A list of strings, where each string is a BED line for a probe.
    """
    chrom, in_file, window_size = args
    probes = []

    try:
        logging.info(f"Processing chromosome: {chrom}")
        # Read the bed file into a pandas DataFrame, filtering for the current chromosome
        df_all = pd.read_csv(
            in_file,
            sep="\t",
            header=None,
            names=["chrom", "start", "end"],
            low_memory=False,
        )

        df = df_all[df_all["chrom"] == chrom].copy()

        # Remove df_all to free memory
        del df_all

        # If there are fewer CpGs than the window size, return empty list
        if df.empty or len(df) < window_size:
            return probes

        # Ensure the data is sorted by coordinate
        df.sort_values(by="start", inplace=True)

        # Loop through the DataFrame in steps of `window_size` to create non-overlapping probes
        for i in range(0, len(df), window_size):
            # The start of the probe is the 'start' coordinate of the current CpG
            start_coord = df.iloc[i]["start"]

            # The end of the probe is the 'end' coordinate of the last CpG in the window
            # Check to make sure the index is not out of bounds
            if i + window_size <= len(df):
                end_coord = df.iloc[i + window_size - 1]["end"]
            else:
                # If the last window is smaller, the end coordinate is the end of the last CpG
                end_coord = df.iloc[len(df) - 1]["end"]

            # Generate the BED line
            probe_name = f"probe_{chrom}_{start_coord}_{end_coord}"
            probes.append(f"{chrom}\t{start_coord}\t{end_coord}\t{probe_name}\n")

    except Exception as e:
        logging.error(f"Error processing chromosome {chrom}: {e}")
        return []

    return probes


def get_chromosomes(file_path):
    """
    Reads the BED file and returns a list of unique chromosome names.
    """
    try:
        with open(file_path, "rt") as f:

            df = pd.read_csv(
                f,
                sep="\t",
                header=None,
                names=["chrom", "start", "end"],
                low_memory=False,
            )
            return df["chrom"].astype(str).unique().tolist()
    except FileNotFoundError:
        logging.error(f"The file '{file_path}' was not found.")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An error occurred while reading chromosomes: {e}")
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        logging.error(
            "Usage: python create_probes.py <input.bed.gz> <output.bed.gz> <window_size> <log_file>"
        )
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    window_size = int(sys.argv[3])
    log = sys.argv[4]

    # Set up logging
    logging.basicConfig(
        format="%(levelname)s:%(asctime)s: %(message)s",
        level=logging.DEBUG,
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.FileHandler(log)],
    )
    logging.info("Starting probe creation process...")
    logging.info(f"Input file: {input_file}")
    logging.info(f"Output file: {output_file}")
    logging.info(f"Window size: {window_size}")

    # Get the list of all unique chromosomes to process
    all_chroms = get_chromosomes(input_file)

    # Set the number of parallel processes
    num_processes = 10

    logging.info(
        f"Starting parallel processing of {len(all_chroms)} chromosomes with a window size of {window_size} CpGs..."
    )

    # Create a pool of worker processes
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Create a list of arguments for each process
        args_list = [(chrom, input_file, window_size) for chrom in all_chroms]

        # Process chromosomes in parallel (no tqdm)
        all_results = pool.map(create_probes, args_list)

    # Sort the results first by chromosome and then by start position
    all_results_flat = [line for result_list in all_results for line in result_list]

    # Extract chromosome and start for sorting
    parsed_results = []
    for line in all_results_flat:
        parts = line.strip().split("\t")
        parsed_results.append((parts[0], int(parts[1]), line))

    parsed_results.sort(key=itemgetter(0, 1))

    # Write the sorted results to the output file
    with open(output_file, "wt") as f:
        for _, _, line in parsed_results:
            f.write(line)

    print(f"Successfully created probes and saved to {output_file}")
