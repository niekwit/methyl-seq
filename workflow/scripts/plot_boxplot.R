# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load required libraries
library(tidyverse)
library(cowplot)
library(data.table)
library(cowplot)
#library(gground)

# Load CpG methylation data
data <- read.delim(snakemake@input[[1]], header = FALSE)
data %>%
  setNames(c("count", "probe", "methylation_status", "region", "condition"))

# Check if there are any probes with only one entry per condition per feature.
# If a probe has only Z or z, the zero values will be missing in data of the
# other methylation call (due to the way data is generated with uniq -c).
# So if a probe has only one entry, add the missing methylation
# call with zero counts.

# Use data.table to complete the data, which is much faster than tidyr::complete

# Convert to data.table
setDT(data)

# Complete the data
completed_data <- dcast(
  data,
  probe + condition + region ~ methylation_status,
  value.var = "count",
  fill = 0
)

# Convert back to long format
completed_data_long <- melt(
  completed_data,
  id.vars = c("probe", "condition", "region"),
  variable.name = "methylation_status",
  value.name = "count"
)

# Calculate percentage methylation per probe
# i.e. percentage of Z of total calls (Z + z)
df <- completed_data_long %>%
  as.data.frame() %>%
  group_by(probe, condition, region) %>%
  mutate(total_calls = sum(count)) %>%
  filter(methylation_status == "Z") %>%
  mutate(methylated_calls = count) %>%
  summarise(
    total_calls = unique(total_calls),
    methylated_calls = unique(methylated_calls),
    perc_methylation = (methylated_calls / total_calls) * 100
  ) %>%
  ungroup()

# Change level order of conditions
sample_info <- read.csv("config/samples.csv", header = TRUE)
condition_levels <- unique(sample_info$condition)

df$condition <- factor(
  df$condition,
  levels = condition_levels
)

# Change level order of regions
region_levels <- snakemake@params[["regions"]]
df$region <- factor(
  df$region,
  levels = region_levels
)

# Create box plots
p <- ggplot(df, aes(x = condition, y = perc_methylation, fill = condition)) +
  geom_boxplot(
    width = 0.5,
    position = position_dodge(0.9),
    outlier.size = 0.5,
  ) +
  facet_wrap(~region) +
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
  width = 6,
  height = 4,
)
