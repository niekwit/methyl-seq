"""
This script aggregates multiple BedGraph files and computes the average score for each genomic interval across all input files. The script is designed
to be used within a Snakemake workflow and expects Snakemake to provide input and output file paths,
as well as a condition wildcard.

Functionality:
- Reads multiple BedGraph files, each containing genomic intervals and associated scores.
- Aggregates scores for each unique (chromosome, start, end) interval across all input files.
- Calculates the average score for each interval, dividing the sum by the total number of input files,
    regardless of whether a particular interval is present in all files (missing values are treated as 0).
- Writes the averaged BedGraph to the specified output file.
- Logs progress and key variables for debugging purposes.

Assumptions:
- BedGraph files are space or tab separated and may contain header/meta lines starting with 'track', 'browser', or '#'.
- Each BedGraph interval is unique per file (e.g., Bismark output: 1-base wide intervals).
- The script is executed within a Snakemake rule, with appropriate input, output, and log file paths.

Inputs (via Snakemake):
- snakemake.input["bg"]: List of input BedGraph file paths.
- snakemake.output["bg"]: Output BedGraph file path.
- snakemake.log[0]: Log file path.
- snakemake.wildcards["condition"]: Condition name for logging.

Output:
- A BedGraph file containing the average score for each genomic interval across all input files.

"""
from collections import defaultdict
import logging

# Set up logging
log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(log)],
)

# Snakemake variables
input_files = snakemake.input["bg"]
condition = snakemake.wildcards["condition"]
output_file = snakemake.output["bg"]

# Log variables for debugging
logging.info(f"Processing condition: {condition}")
logging.info(f"Input bedGraph files: {input_files}")
logging.info(f"Condition: {condition}")
logging.info(f"Output bedGraph files: {output_file}")


### Use a dictionary of dictionaries to store scores:
# { (chrom, start, end): [score_sample1, score_sample2, ...] }
all_scores = defaultdict(lambda: [0.0] * len(input_files))

logging.info(f"Reading and aggregating data from {len(input_files)} files...")

### Aggregate Scores
for i, file_path in enumerate(input_files):
    logging.info(f"Processing {file_path}...")

    # Note: BedGraph files are typically space or tab separated
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith(("track", "browser", "#")):
                continue  # Skip header/meta lines

            parts = line.strip().split()
            if len(parts) >= 4:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                score = float(parts[3])

                # Bismark BedGraphs are 1-base wide, so start/end are unique
                key = (chrom, start, end)

                # Store the score in the correct list position (index i)
                all_scores[key][i] = score

logging.info("Aggregation complete. Calculating means...")

### Calculate Mean and Write Output
output_lines = []
for (chrom, start, end), scores in all_scores.items():

    # Calculate the total sum of all scores (missing samples contribute 0.0)
    sum_scores = sum(scores)

    # Check if the site was covered in AT LEAST ONE sample (i.e., the total score > 0)
    if sum_scores > 0:

        # Divide the sum by the TOTAL number of input files (3 in your case)
        # This achieves the 100/3 = 33.3333 calculation
        total_replicates = len(scores)
        average_score = sum_scores / total_replicates

        # Format the output line
        output_lines.append(f"{chrom}\t{start}\t{end}\t{average_score:.2f}")

# Write all results to the output file
with open(output_file, "w") as out_f:
    out_f.write("\n".join(output_lines) + "\n")
logging.info(f"Successfully wrote average BedGraph to {output_file}")
