# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

te_names <- snakemake@params[["te_names"]]
reference_condition <- snakemake@params[["reference_condition"]]
ko_conditions <- snakemake@params[["ko_conditions"]]
cutoffs <- setNames(snakemake@params[["cutoffs"]], te_names)

out_dir <- "results/te_5utr_status"
bed_col_names <- c("chrom", "start", "end", "name", "score", "strand")

write_empty_bed <- function(te, ko_condition, element_status) {
  write_tsv(
    data.frame(matrix(ncol = 6, nrow = 0)),
    file.path(out_dir, paste0(te, "_", ko_condition, "_", element_status, ".bed")),
    col_names = FALSE
  )
}

# bigWigAverageOverBed output (name, size, covered, sum, mean0, mean), one
# row per element's 5' UTR x TE x condition -- see line1_5utr_bigwig_average
# in te_utr.smk. "name" is the element's own "_5UTR" feature id (see
# line1_annotate_utrs_orfs.py); strip that suffix to recover the element's
# base id, shared with its "_remainder" row in bed/{te_name}_5UTR_remainder.bed.
data <- read.delim(
  snakemake@input[["avg"]],
  header = FALSE,
  col.names = c(
    "name", "size", "covered", "sum", "mean0", "mean", "TE", "condition"
  )
) %>%
  mutate(base_id = sub("_5UTR$", "", name))

if (nrow(data) == 0) {
  # No near-full-length L1 elements were confidently annotated for any
  # configured TE -- data-dependent (sequence divergence/completeness),
  # not necessarily a config error, so render a placeholder page and
  # empty BED files instead of erroring.
  message(
    "No 5' UTR methylation data -- writing placeholder plot and empty BED files"
  )
  p <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No near-full-length L1 elements were confidently annotated"
    ) +
    theme_void()
  ggsave(filename = snakemake@output[["pdf"]], plot = p, width = 6, height = 4)
  write.csv(data, file = snakemake@output[["csv"]], row.names = FALSE)

  for (te in te_names) {
    for (ko in ko_conditions) {
      write_empty_bed(te, ko, "active")
      write_empty_bed(te, ko, "inactive")
    }
  }

  sink(type = "message")
  sink(type = "output")
  quit(save = "no", status = 0)
}

# Condition order from config/samples.csv, matching plot_te_5utr_boxplot.R
sample_info <- read.csv("config/samples.csv", header = TRUE)
condition_levels <- unique(sample_info$condition)
data$condition <- factor(data$condition, levels = condition_levels)
data$TE <- factor(data$TE, levels = te_names)

# Histogram of mean 5' UTR methylation per element, one facet per TE,
# coloured by condition (reference + every KO condition), with each TE's
# configured cutoff drawn as a reference line
cutoff_df <- data.frame(
  TE = factor(te_names, levels = te_names),
  cutoff = unname(cutoffs[te_names])
)

p <- ggplot(data, aes(x = mean, fill = condition)) +
  geom_histogram(
    position = "identity",
    alpha = 0.6,
    bins = 50,
    colour = NA
  ) +
  geom_vline(
    data = cutoff_df,
    aes(xintercept = cutoff),
    linetype = "dashed",
    colour = "black"
  ) +
  facet_wrap(vars(TE), scales = "free_y") +
  labs(
    x = "Mean 5' UTR CpG methylation (%)",
    y = "Count",
    fill = NULL
  ) +
  theme_cowplot(14) +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "top"
  )

ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = 3 + length(te_names) * 2.5,
  height = 5
)

write.csv(data, file = snakemake@output[["csv"]], row.names = FALSE)

# --- Classify each element's 5' UTR methylation status, per TE, for each
# non-reference ("KO") condition against the single reference condition
# (config DMR:reference_condition):
#   active:   reference >= cutoff & KO < cutoff (lost methylation upon KO)
#   inactive: reference >= cutoff & KO >= cutoff (stays methylated in KO)
# Elements with reference < cutoff are excluded from both (never started
# out methylated, so "loss upon KO" doesn't apply). Each classified
# element's full genomic span (5' UTR + remainder rows collapsed into one)
# is written out as an active/inactive BED file -- one pair per TE x KO
# condition, declared as this rule's active/inactive outputs.
for (i in seq_along(te_names)) {
  te <- te_names[i]
  # snakemake@input[["beds"]] is expand()'d over te_names in the same
  # order (see plot_line1_5utr_methylation_status in te_utr.smk)
  bed <- read.delim(
    snakemake@input[["beds"]][i],
    header = FALSE,
    col.names = bed_col_names
  ) %>%
    mutate(base_id = sub("_(5UTR|remainder)$", "", name))

  cutoff <- cutoffs[[te]]

  # One row per element, one column per condition's mean 5' UTR methylation
  wide <- data %>%
    filter(TE == te) %>%
    select(base_id, condition, mean) %>%
    pivot_wider(names_from = condition, values_from = mean)

  for (ko in ko_conditions) {
    status <- wide %>%
      filter(
        !is.na(.data[[reference_condition]]),
        !is.na(.data[[ko]])
      ) %>%
      transmute(
        base_id,
        status = case_when(
          .data[[reference_condition]] >= cutoff & .data[[ko]] < cutoff ~
            "active",
          .data[[reference_condition]] >= cutoff & .data[[ko]] >= cutoff ~
            "inactive",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(status))

    for (element_status in c("active", "inactive")) {
      ids <- status %>%
        filter(status == element_status) %>%
        pull(base_id)

      if (length(ids) == 0) {
        write_empty_bed(te, ko, element_status)
        next
      }

      subset_bed <- bed %>%
        filter(base_id %in% ids) %>%
        group_by(chrom, base_id, strand) %>%
        summarise(start = min(start), end = max(end), .groups = "drop") %>%
        transmute(chrom, start, end, name = base_id, score = ".", strand) %>%
        arrange(chrom, start)

      write_tsv(
        subset_bed,
        file.path(out_dir, paste0(te, "_", ko, "_", element_status, ".bed")),
        col_names = FALSE
      )
    }
  }
}

sink(type = "message")
sink(type = "output")
