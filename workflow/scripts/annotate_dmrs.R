# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Mm.eg.db)
library(AnnotationDbi)

# Get input files from Snakemake
dmr_hyper_file <- snakemake@input[["hyper"]]
dmr_hypo_file <- snakemake@input[["hypo"]]
diff_tiles_file <- snakemake@input[["diff_tiles"]]

# Get output files from Snakemake
dmr_hyper_annotated_output <- snakemake@output[["hyper"]]
dmr_hypo_annotated_output <- snakemake@output[["hypo"]]
distance_plot <- snakemake@output[["distance"]]
distribution_plot <- snakemake@output[["distribution"]]
volcano_plot <- snakemake@output[["volcano"]]

# Get Snakemake parameters
genome <- snakemake@config[["genome"]]
difference_threshold <- snakemake@config[["difference_threshold"]]
qvalue_threshold <- snakemake@config[["qvalue_threshold"]]

# Load bed files and convert to GRanges
hyper_gr <- rtracklayer::import(dmr_hyper_file, format = "bed")
hypo_gr <- rtracklayer::import(dmr_hypo_file, format = "bed")

# Define the TxDb
if (genome == "mm39") {
  txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
} else if (genome == "hg38") {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
} else {
  stop("Unsupported genome specified in config.yaml")
}

# Convert to Ensembl style
seqlevelsStyle(txdb) <- "Ensembl"

# Annotate DMRs using ChIPseeker
gr_list <- list(
  Hypermethylated = hyper_gr,
  Hypomethylated = hypo_gr
)
peakAnnoList <- lapply(
  gr_list,
  annotatePeak,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  verbose = FALSE
)

# Map Entrez IDs to Symbols for both Hyper and Hypo objects
peakAnnoList <- lapply(peakAnnoList, function(anno_obj) {
  # Extract the annotation metadata
  anno_df <- as.data.frame(anno_obj@anno)

  # Query the database for Symbols and Ensembl IDs
  gene_map <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(anno_df$geneId),
    columns = c("SYMBOL", "ENSEMBL", "GENENAME"),
    keytype = "ENTREZID"
  )

  # Merge the names back into the annotation data
  # Use dplyr for a clean join
  updated_anno <- left_join(anno_df, gene_map, by = c("geneId" = "ENTREZID"))

  # Convert back to GRanges and re-insert into the object
  # We must keep the same metadata structure
  updated_gr <- makeGRangesFromDataFrame(
    updated_anno,
    keep.extra.columns = TRUE
  )

  # Re-assign the seqinfo to keep the 'mm39' genome metadata
  seqinfo(updated_gr) <- seqinfo(anno_obj@anno)
  anno_obj@anno <- updated_gr

  return(anno_obj)
})

# Add column to mark if genes are unique to hyper or hypo group
peakAnnoList_df <- lapply(peakAnnoList, as.data.frame)
for (i in seq_along(peakAnnoList_df)) {
  tmp <- peakAnnoList_df[[i]]

  # Mark genes as "Unique to Hyper" or "Unique to Hypo"
  label <- names(peakAnnoList_df)[i]
  if (label == "Hypermethylated") {
    tmp$Group <- ifelse(
      tmp$geneId %in% genes[["Hypermethylated"]],
      "Unique to Hyper",
      "Shared"
    )
  } else if (label == "Hypomethylated") {
    tmp$Group <- ifelse(
      tmp$geneId %in% genes[["Hypomethylated"]],
      "Unique to Hypo",
      "Shared"
    )
  }
  peakAnnoList_df[[i]] <- tmp
}

# Save annotated peaks as csv
write.csv(
  peakAnnoList_df$Hypermethylated,
  dmr_hyper_annotated_output,
  row.names = FALSE
)
write.csv(
  peakAnnoList_df$Hypomethylated,
  dmr_hypo_annotated_output,
  row.names = FALSE
)

# Create distribution plot
pdf(distribution_plot, width = 6, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()

# Create distance to TSS plot
pdf(distance_plot, width = 6, height = 3)
plotDistToTSS(peakAnnoList)
dev.off()

# Create volcano plot
diff_tiles <- readRDS(diff_tiles_file)

diff_tiles_df <- as.data.frame(diff_tiles_noLINE1_gr) %>%
  mutate(
    sig = case_when(
      qvalue < qvalue_threshold &
        meth.diff > difference_threshold ~ "Hypermethylated",
      qvalue < qvalue_threshold &
        meth.diff < -difference_threshold ~ "Hypomethylated",
      TRUE ~ "Not significant"
    )
  )

p_volcano <- ggplot(
  diff_tiles_df,
  aes(x = meth.diff, y = -log10(qvalue), colour = sig)
) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(
    values = c(
      "Hypermethylated" = "#E41A1C",
      "Hypomethylated" = "#377EB8",
      "Not significant" = "grey70"
    )
  ) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "black") +
  geom_vline(xintercept = c(-25, 25), linetype = "dashed", colour = "black") +
  labs(
    x = "Methylation difference (%)",
    y = expression(-log[10](q - value)),
    colour = NULL
  ) +
  theme_cowplot(12) +
  theme(legend.position = "bottom")

ggsave(volcano_plot, p_volcano, width = 6, height = 4, bg = "white")
