# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")


library(tidyverse)
library(cowplot)

preseq_files <- snakemake@input[["preseq"]]
pdf_files <- snakemake@output[["pdf"]]

# df to store efficiency data and recommendation logic
df <- data.frame(
  sample = character(),
  efficiency = numeric(),
  recommendation = character(),
  stringsAsFactors = FALSE
)

# Process and plot preseq data
for (file in preseq_files) {
  # Get sample name from file path
  sample <- gsub(".txt$", "", basename(file))

  # Read data
  preseq_data <- read.delim(file, header = TRUE) %>%
    filter(if_all(everything(), ~ !is.na(.)))

  # Calculate the discovery efficiency at the maximum projected point
  latest_point <- preseq_data %>% filter(TOTAL_READS == max(TOTAL_READS))
  efficiency <- latest_point$EXPECTED_DISTINCT / latest_point$TOTAL_READS

  # Define recommendation logic
  recommendation <- case_when(
    efficiency > 0.8 ~ "High Complexity: Recommend deeper sequencing.",
    efficiency > 0.5 ~ "Moderate Saturation: Diminishing returns starting.",
    TRUE ~ "Saturated: Sequencing deeper is not recommended."
  )

  status_color <- case_when(
    efficiency > 0.8 ~ "#28a745", # Green
    efficiency > 0.5 ~ "#ffc107", # Amber
    TRUE ~ "#dc3545" # Red
  )

  # Create the plot
  p <- ggplot(preseq_data, aes(x = TOTAL_READS, y = EXPECTED_DISTINCT)) +
    # Add confidence interval ribbon
    geom_ribbon(
      aes(ymin = LOWER_0.95CI, ymax = UPPER_0.95CI),
      fill = "grey80",
      alpha = 0.5
    ) +
    # Add the 'Ideal' line (y = x, where every read is unique)
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "red",
      alpha = 0.5
    ) +
    # Add the complexity curve
    geom_line(color = "#007bff", size = 1.2) +
    geom_point(color = "#007bff", size = 2) +
    scale_x_continuous(
      labels = label_number(scale = 1e-6, suffix = "M"),
      guide = guide_axis(cap = "both")
    ) +
    scale_y_continuous(
      labels = label_number(scale = 1e-6, suffix = "M"),
      guide = guide_axis(cap = "both"),
      limits = c(0, max(preseq_data$EXPECTED_DISTINCT) * 1.2),
    ) +
    labs(
      title = "Preseq Library Complexity Projection",
      subtitle = paste0(
        "Efficiency at max projection: ",
        round(efficiency * 100, 1),
        "%"
      ),
      x = "Total Reads Sequenced (Millions)",
      y = "Expected Distinct Fragments (Millions)"
    ) +
    theme_cowplot(10) +
    # Visual Recommendation Box
    annotate(
      "label",
      x = max(preseq_data$TOTAL_READS) * 0.05,
      y = max(preseq_data$EXPECTED_DISTINCT) * 0.9,
      label = paste("Recommendation:\n", recommendation),
      fill = status_color,
      color = "white",
      fontface = "bold",
      hjust = 0,
      size = 2
    )

  # Get the corresponding PDF output file path
  pdf_file <- pdf_files[grep(sample, pdf_files)]

  # Save the plot
  ggsave(pdf_file, p, height = 5, width = 6)

  # Append to summary dataframe
  df <- rbind(
    df,
    data.frame(
      sample = sample,
      efficiency = round(efficiency, 3),
      recommendation = recommendation,
      stringsAsFactors = FALSE
    )
  )
}
