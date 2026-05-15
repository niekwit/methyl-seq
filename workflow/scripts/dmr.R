# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(data.table)
library(methylKit)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

# Get input files from Snakemake
cpgot_files <- snakemake@input[["cpgot"]]
cpgob_files <- snakemake@input[["cpgob"]]

# Get paramaters from Snakemake
samples <- snakemake@params[["samples"]]
reference_condition <- snakemake@params[["ref_cond"]]
genome <- snakemake@config[["genome"]]
tile_size <- snakemake@params[["tile_size"]]
step_size <- snakemake@params[["step_size"]]
min_per_group <- strtoi(snakemake@params[["min_per_group"]])
difference_threshold <- as.numeric(snakemake@params[["difference_threshold"]])
qvalue_threshold <- as.numeric(snakemake@params[["qvalue_threshold"]])

# Get output files from Snakemake
dmr_hyper_output <- snakemake@output[["hyper"]]
dmr_hypo_output <- snakemake@output[["hypo"]]
dmr_tiles_output <- snakemake@output[["meth_tiles"]]
dmr_diff_tiles_output <- snakemake@output[["diff_tiles"]]
dmr_diff_tiles_sig_output <- snakemake@output[["sign_diff_tiles"]]


process_bismark_data <- function(
  ot_file,
  ob_file,
  assembly
) {
  sample_id <- str_replace(basename(ot_file), "_OT_CpG.txt", "")
  message(paste("Processing sample:", sample_id))

  # Read only the necessary columns (Chr, Start, Call)
  cols <- c("V3", "V4", "V5")

  ot <- fread(ot_file, select = cols, col.names = c("chr", "start", "call"))
  ob <- fread(ob_file, select = cols, col.names = c("chr", "start", "call"))

  # Combine strands into one data.table
  dt <- rbindlist(list(ot, ob))
  rm(ot, ob)
  gc() # Garbage collection to keep memory lean

  # Collapse i and i+1 into a single CpG unit
  # We round even positions down to the preceding odd position.
  # This aligns the C (forward) and G (reverse) to the same coordinate.
  dt[, start := ifelse(start %% 2 == 0, start - 1, start)]

  # Aggregate by position to calculate coverage and methylation counts
  summary_dt <- dt[,
    .(
      coverage = .N,
      numCs = sum(call == "Z"), # Z = Methylated
      numTs = sum(call == "z") # z = Unmethylated
    ),
    by = .(chr, start)
  ]

  # Format for methylKit requirements
  summary_dt[, `:=`(
    end = start,
    strand = "*"
  )]

  # Ensure column order matches methylKit methylRaw standard
  setcolorder(
    summary_dt,
    c("chr", "start", "end", "strand", "coverage", "numCs", "numTs")
  )

  # Return as methylRaw object
  return(new(
    "methylRaw",
    as.data.frame(summary_dt),
    sample.id = sample_id,
    assembly = assembly,
    context = "CpG",
    resolution = "base"
  ))
}

# Process each sample's OT and OB files into methylRaw objects
meth_list <- mapply(
  process_bismark_data,
  ot_file = cpgot_files,
  ob_file = cpgob_files,
  MoreArgs = list(assembly = genome),
  SIMPLIFY = FALSE
)

# Create a vector that contains 1 or 0 for each sample,
# indicating whether it belongs to the treatment (1) condition or not
ref_vector <- as.integer(!str_detect(samples, reference_condition))

# Create the methylRawList
methylobj <- new("methylRawList", meth_list, treatment = ref_vector)
meth_final <- methylKit::unite(
  methylobj,
  destrand = FALSE,
  min.per.group = min_per_group
)

# Tiling windows for DMR analysis
tiles <- tileMethylCounts(
  methylobj,
  win.size = tile_size,
  step.size = step_size
)
meth_tiles <- methylKit::unite(tiles, min.per.group = min_per_group)

# Save all tiles to file
saveRDS(meth_tiles, file = dmr_tiles_output)

# Get differentially methylated tiles
diff_tiles <- calculateDiffMeth(meth_tiles, mc.cores = 12)
saveRDS(diff_tiles, file = dmr_diff_tiles_output)

# Convert differential tiles to GRanges
# Filter for significant hits first based on thresholds
diff_tiles_sig <- getMethylDiff(
  diff_tiles,
  difference = difference_threshold,
  qvalue = qvalue_threshold
)
saveRDS(diff_tiles_sig, file = dmr_diff_tiles_sig_output)

tiles_gr <- as(diff_tiles_sig, "GRanges")

# Split into hyper- and hypo-methylated DMRs and save as BED format
hyper_gr <- tiles_gr[tiles_gr$meth.diff > 0]
rtracklayer::export(hyper_gr, dmr_hyper_output)

hypo_gr <- tiles_gr[tiles_gr$meth.diff < 0]
rtracklayer::export(hypo_gr, dmr_hypo_output)
