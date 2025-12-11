# Nurefsan Sariipek and Peter van Galen, 250713
# Visualize the Numbat results files (see https://kharchenkolab.github.io/numbat/articles/results.html)

library(tidyverse)
library(numbat)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/09_Numbat"))

# Clear environment variables
rm(list = ls())

# Map IDs. This is because Numbat was run before we changed all other analyses to a new ID
# fmt: skip
id_map <- tribble(
  ~OLD_ID, ~NEW_ID,
  "P05", "P20",
  "P07", "P22",
  "P08", "P23",
  "P09", "P30",
  "P10", "P31",
  "P12", "P33"
) %>% column_to_rownames("NEW_ID")

# Loop to plot copy number landscape
for (pt in rownames(id_map)) {
  #pt <- "P20"
  print(paste("Analyzing", pt))

  # Load Numbat object (need to mount Broad storage first)
  nb <- Numbat$new(
    out_dir = paste0(
      "/Volumes/broad_vangalenlab/sariipek/numbat/",
      id_map[pt, "OLD_ID"]
    )
  )

  # Generate copy number landscape plot
  cnv_plot <- nb$plot_phylo_heatmap(
    clone_bar = TRUE,
    p_min = 0.9,
    pal_clone = c(
      "1" = "gray",
      "2" = "#377EB8",
      "3" = "#4DAF4A",
      "4" = "#984EA3"
    )
  )

  #cnv_plot
  # Aggregate cells by clone and visualize CNV events in pseudobulks
  bulk_plot <- nb$bulk_clones %>%
    filter(n_cells > 50) %>%
    plot_bulks(legend = TRUE)

  # Save as pdf
  pdf(paste0("9.3_", pt, "_results.pdf"), width = 12, height = 4)
  print(cnv_plot)
  print(bulk_plot)
  dev.off()
}
