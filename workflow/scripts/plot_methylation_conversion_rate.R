# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load required libraries
library(tidyverse)
library(cowplot)
library(data.table)

# Load coverage data
cov_files <- snakemake@input[["cov"]]
cov_list <- lapply(cov_files, fread, header = FALSE)
cov_data <- rbindlist(cov_list)
colnames(cov_data) <- c(
  "count",
  "contig",
  "pos",
  "methylation_status",
  "sample"
)

# Check if there are any probes with only one entry per condition per feature.
# If a probe has only Z or z, the zero values will be missing in data of the
# other methylation call (due to the way data is generated with uniq -c).
# So if a probe has only one entry, add the missing methylation
# call with zero counts.

# Use data.table to complete the data, which is much faster than tidyr::complete

# Convert to data.table
setDT(cov_data)

# Complete the data
completed_data <- dcast(
  cov_data,
  contig + sample + pos ~ methylation_status,
  value.var = "count",
  fill = 0
)

# For each contig and sample, sum the counts of Z and z
completed_data <- completed_data %>%
  as.data.frame() %>%
  group_by(contig, sample) %>%
  mutate(
    total_Z = sum(Z, na.rm = TRUE),
    total_z = sum(z, na.rm = TRUE),
    methylation_rate = total_Z / (total_Z + total_z) * 100
  ) %>%
  ungroup()

# Load sample info to set factor levels
sample_info <- read.csv("config/samples.csv")
completed_data$sample <- factor(
  completed_data$sample,
  levels = sample_info$sample
)

# Plot methylation conversion rate
# Plot bar side-by-side for different contigs
p <- ggplot(
  completed_data,
  aes(x = sample, y = methylation_rate, fill = contig)
) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    x = NULL,
    y = "Methylation Conversion Rate (%)"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme_cowplot() +
  scale_fill_manual(
    values = c("phage_lambda" = "#1f77b4", "plasmid_puc19c" = "#ff7f0e")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot
ggsave(snakemake@output[["pdf"]], plot = p)

# Save data
write.csv(
  completed_data,
  snakemake@output[["csv"]],
  row.names = FALSE,
  quote = FALSE
)
