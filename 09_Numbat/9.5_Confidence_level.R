# Nurefsan Sariipek, 250402
# Check the confidence levels of Numbat calls

# Load the libraries
library(tidyverse)
library(Seurat)

# Empty environment
rm(list = ls())

# Set working directory (for Nurefsan)
setwd("~/TP53_ImmuneEscape/09_Numbat/")

# For Peter
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/09_Numbat")

# Favorite function
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Load the saved Seurat object
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Get the Numbat clone files
tsv_files <- list.files("Numbat_calls/", pattern = "*.tsv", full.names = TRUE)

# Use read_delim with tab delimiter
all_pcnv <- tsv_files %>%
  map_dfr(
    ~ read_delim(.x, delim = "\t", show_col_types = FALSE) %>%
      select(cell, p_cnv_x, p_cnv_y, p_cnv, compartment_opt) %>%
      mutate(sample = tools::file_path_sans_ext(basename(.x)))
  )

# Add paient ID
all_pcnv$patient_id <- cutf(all_pcnv$sample, d = "_")

# Check result
head(all_pcnv)

# Density plots for confidence
# p_cnv_x
ggplot(all_pcnv, aes(x = p_cnv_x, fill = patient_id)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Posterior probability per sample (expression)",
    x = "Posterior probability",
    y = "Density"
  ) +
  facet_wrap(~patient_id, scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("8.5.1_Expression_pcnv.pdf", width = 12, height = 8)

# p_cnv_y
ggplot(all_pcnv, aes(x = p_cnv_y, fill = patient_id)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Posterior probability per sample (allele)",
    x = "Posterior probability",
    y = "Density"
  ) +
  facet_wrap(~patient_id, scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("8.5.2_Allele_pcnv.pdf", width = 12, height = 8)

# p_cnv
ggplot(all_pcnv, aes(x = p_cnv, fill = patient_id)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Posterior probability per sample (p_cnv)",
    x = "Posterior probability",
    y = "Density"
  ) +
  facet_wrap(~patient_id, scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("8.5.3_Posterior_pcnv.pdf", width = 12, height = 8)

# To merge with Seurat data, use common cell barcodes
all_pcnv$barcode <- paste0(all_pcnv$patient_id, "_", all_pcnv$cell)
seu$new_barcode <- paste0(
  seu$patient_id,
  "_",
  cutf(colnames(seu), d = "_", f = 3)
)
all(all_pcnv$barcode %in% seu$new_barcode)

# Subset & merge
seu_subset <- subset(seu, new_barcode %in% all_pcnv$barcode)
seu_subset <- RenameCells(seu_subset, new.names = seu_subset$new_barcode)
seu_subset <- AddMetaData(
  seu_subset,
  column_to_rownames(all_pcnv, var = "barcode")
)

# Plot the UMAP with probability
FeaturePlot(seu_subset, features = "p_cnv_y") +
  scale_color_gradient(low = "lightgrey", high = "firebrick") +
  labs(
    title = "Posterior CNV probability (p_cnv_y)",
    color = "CNV\nProbability"
  ) +
  theme(
    aspect.ratio = 1,
    axis.line = element_blank(),
    panel.border = element_rect(color = "black")
  )

ggsave("8.5.4_UMAP_cnv_prob.pdf", width = 7, height = 7)


# Show that T cells that were identified as tumor cells have low confidence level
seu_subset2 <- subset(seu_subset, numbat_compartment == "tumor")
seu_subset2$cnv_confidence <- ifelse(
  seu_subset2$p_cnv > 0.9,
  yes = "High",
  no = "Low"
)

DimPlot(seu_subset2, group.by = "cnv_confidence", pt.size = 0.5, alpha = 0.8) +
  scale_color_manual(values = c("Low" = "grey", "High" = "darkred")) +
  labs(
    title = "Tumor cells by CNV call confidence",
    color = "Confidence>0.9"
  ) +
  theme(
    aspect.ratio = 1,
    axis.line = element_blank(),
    panel.border = element_rect(color = "black")
  )

ggsave("8.5.5_Tumor_cells_CNV_confidence.pdf", width = 7, height = 7)
