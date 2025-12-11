# Nurefsan Sariipek and Peter van Galen, 250722
# Pie chart for Figure 1

# Load libraries
library(tidyverse)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/10_Miscellaneous"))

# Clear environment variables
rm(list = ls())

# Define cohort colors
cohort_colors <- c(
  "Remission_MUT" = "#546fb5",
  "Remission_WT" = alpha("#546fb5", 0.4),
  "Relapse_WT" = alpha("#e54c35", 0.4),
  "Relapse_MUT" = "#e54c35"
)

# Create the data
df <- tibble(
  cohort = c("Remission_MUT", "Remission_WT", "Relapse_WT", "Relapse_MUT"),
  n = c(6, 13, 6, 4)
) %>%
  mutate(cohort = factor(cohort, levels = cohort)) %>%
  mutate(
    label = paste0(cohort, "\n(n=", n, ")")
  )

# Calculate positions for text
df <- df %>%
  mutate(csum = cumsum(n), pos = csum - n / 2) %>%
  mutate(pos = sum(n) - pos)

# Plot
pie_chart <- df %>%
  ggplot(aes(x = "", y = n, fill = cohort)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(y = pos, label = label), color = "black", size = 3) +
  scale_fill_manual(values = cohort_colors) +
  theme_void() +
  theme(legend.position = "none")

# View
pie_chart

# Save
pdf("10.1_Piechart.pdf", width = 5, height = 5)
pie_chart
dev.off()
