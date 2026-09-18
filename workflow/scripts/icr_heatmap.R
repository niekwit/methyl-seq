# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

# bigWigAverageOverBed output (name, size, covered, sum, mean0, mean), one
# row per ICR per condition, concatenated across conditions with a
# "condition" column appended -- see the icr_scores rule in icr.smk
df <- read.delim(
  snakemake@input[["scores"]],
  header = FALSE,
  col.names = c("name", "size", "covered", "sum", "mean0", "mean", "condition")
)

# Region order matches the input BED file's own row order (chrom, start,
# end, name)
icr_bed <- read.delim(snakemake@input[["icr_bed"]], header = FALSE)
name_levels <- icr_bed[[4]]
df$name <- factor(df$name, levels = name_levels)

# Condition order from config/samples.csv, matching plot_boxplot.R
sample_info <- read.csv("config/samples.csv", header = TRUE)
condition_levels <- unique(sample_info$condition)
df$condition <- factor(df$condition, levels = condition_levels)

text_colour <- "#ffffff"

# Heatmap of average %CpG methylation per ICR per condition
p <- ggplot(df, aes(x = condition, y = name, fill = mean)) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "F",
    direction = -1,
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20)
  ) +
  labs(
    x = NULL,
    y = "Paternally imprinted\nregions (ICRs)",
    fill = "Average %CpG\nmethylation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.y = element_text(face = "italic")
  ) +
  geom_text(aes(label = round(mean, 1)), color = text_colour, size = 4)

ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = 1.5 + length(condition_levels) * 0.9,
  height = 1.5 + length(name_levels) * 0.5
)

write.csv(df, file = snakemake@output[["csv"]], row.names = FALSE)
