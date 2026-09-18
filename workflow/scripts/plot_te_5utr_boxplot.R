# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load required libraries
library(tidyverse)
library(cowplot)
library(data.table)

te_names <- snakemake@params[["te_names"]]
min_reads <- snakemake@params[["min_reads"]]

# Load CpG methylation call counts (count, region_id, methylation_status,
# TE, condition), one row per element-UTR/remainder x methylation-call-type
# x TE x condition -- see line1_utr_boxplot_data in te_utr.smk
data <- read.delim(snakemake@input[[1]], header = FALSE) %>%
  setNames(c("count", "region_id", "methylation_status", "TE", "condition")) %>%
  # 5' UTR vs. remainder is encoded in the feature name's suffix (see
  # line1_annotate_utrs_orfs.py / generate_line1_utr_regions_* rules)
  mutate(region = ifelse(grepl("_5UTR$", region_id), "5' UTR", "L1 remainder"))

if (nrow(data) == 0) {
  # No near-full-length L1 elements were confidently annotated for any
  # configured TE -- data-dependent (sequence divergence/completeness),
  # not necessarily a config error, so render a placeholder page instead
  # of erroring.
  message("No 5'UTR/remainder methylation data -- writing placeholder plot")
  p <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No near-full-length L1 elements were confidently annotated"
    ) +
    theme_void()
  ggsave(filename = snakemake@output[["pdf"]], plot = p, width = 6, height = 4)
  write.csv(
    data.frame(
      region_id = character(),
      condition = character(),
      TE = character(),
      region = character(),
      total_calls = integer(),
      methylated_calls = integer(),
      perc_methylation = numeric()
    ),
    file = snakemake@output[["csv"]],
    row.names = FALSE
  )
  sink(type = "message")
  sink(type = "output")
  quit(save = "no", status = 0)
}

# Check if there are any elements with only one entry per condition per
# feature. If an element has only Z or z, the zero values will be missing
# in data of the other methylation call (due to the way data is generated
# with uniq -c). So if an element has only one entry, add the missing
# methylation call with zero counts.

# Use data.table to complete the data, which is much faster than tidyr::complete
setDT(data)

completed_data <- dcast(
  data,
  region_id + condition + TE + region ~ methylation_status,
  value.var = "count",
  fill = 0
)

completed_data_long <- melt(
  completed_data,
  id.vars = c("region_id", "condition", "TE", "region"),
  variable.name = "methylation_status",
  value.name = "count"
)

# Calculate percentage methylation per element
# i.e. percentage of Z of total calls (Z + z)
df <- completed_data_long %>%
  as.data.frame() %>%
  group_by(region_id, condition, TE, region) %>%
  mutate(total_calls = sum(count)) %>%
  filter(methylation_status == "Z") %>%
  mutate(methylated_calls = count) %>%
  summarise(
    total_calls = unique(total_calls),
    methylated_calls = unique(methylated_calls),
    perc_methylation = (methylated_calls / total_calls) * 100
  ) %>%
  ungroup() %>%
  filter(total_calls >= min_reads)

# Change level order of conditions
sample_info <- read.csv("config/samples.csv", header = TRUE)
condition_levels <- unique(sample_info$condition)
df$condition <- factor(df$condition, levels = condition_levels)

# Change level order of TE / region
df$TE <- factor(df$TE, levels = te_names)
df$region <- factor(df$region, levels = c("5' UTR", "L1 remainder"))

# Create box plots: %CpG methylation in the 5' UTR vs. the remainder,
# faceted by TE (family total + configured subfamilies)
p <- ggplot(df, aes(x = condition, y = perc_methylation, fill = condition)) +
  geom_boxplot(
    width = 0.5,
    position = position_dodge(0.9),
    outlier.size = 0.5,
  ) +
  facet_grid(cols = vars(TE), rows = vars(region)) +
  labs(
    x = NULL,
    y = "CpG methylation (%)"
  ) +
  theme_cowplot() +
  theme(
    panel.border = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
  ) +
  scale_fill_manual(values = c("#cccccc", "#5373c6"))

# Save plot
ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = 2 + length(te_names) * 1.5,
  height = 6,
)

# Save data to csv
write.csv(df, file = snakemake@output[["csv"]], row.names = FALSE)
