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

# Define the correct column names for your custom 5-column format
COV_COL_NAMES <- c("count", "contig", "pos", "methylation_status", "sample")

# Read the data using lapply and fread
cov_list <- lapply(
  cov_files,
  data.table::fread,
  header = FALSE,
  col.names = COV_COL_NAMES
)

# Combine the list into one data table
cov_data <- data.table::rbindlist(cov_list)

# Calculate methylation conversion rate per sample and contig
df <- cov_data %>%
  group_by(contig, sample, methylation_status) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = methylation_status,
    values_from = total,
    values_fill = 0
  ) %>%
  mutate(
    total_calls = Z + z,
    methylation_rate = (Z / total_calls) * 100
  ) %>%
  ungroup()

# Load sample info to set factor levels
sample_info <- read.csv("config/samples.csv")
df$sample <- factor(
  df$sample,
  levels = sample_info$sample
)


# Plot methylation conversion rate
p <- ggplot(
  df,
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
  df,
  snakemake@output[["csv"]],
  row.names = FALSE,
  quote = FALSE
)
