# Nurefsan Sariipek, 250402

#Check the confidence levels of Numbat calls
# Loaad the libraries
library(tidyverse)
library(ggplot2)
library(Seurat)

# Empty environment
rm(list=ls())

# Set working directory
setwd("~/TP53_ImmuneEscape/8_Numbat/")
# For Peter (local)
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/8_Numbat")

# Load the saved Seurat objects
seu_combined <- readRDS("~/250505_numbat_combined_seurat.rds")

# Get the numbat clone files
tsv_files <- list.files("Numbat_Calls/", pattern = "*.tsv", full.names = TRUE)

# Use read_delim with tab delimiter
all_pcnv <- tsv_files %>%
  map_dfr(~ read_delim(.x, delim = "\t", show_col_types = FALSE) %>%
            select(cell, p_cnv_x, p_cnv_y, p_cnv, compartment_opt) %>%
            mutate(sample = tools::file_path_sans_ext(basename(.x))))

# Wrangle the data
all_pcnv <- all_pcnv %>%
  mutate(patient_id = str_extract(sample, "^P\\d+"))

# Check result
head(all_pcnv)

# Density plots
# p_cnv_x
ggplot(all_pcnv, aes(x = p_cnv_x, fill = patient_id)) +
  geom_density(alpha = 0.5) +
  labs(title = "Posterior Probability per Sample (expression)",
       x = "Posterior Probability",
       y = "Density") +
  facet_wrap(~patient_id, scales = "free_y") +
  theme_minimal()
ggsave("8.6_expression_pcnv.png", width = 12, height = 8)

# p_cnv_y
ggplot(all_pcnv, aes(x = p_cnv_y, fill = patient_id)) +
  geom_density(alpha = 0.5) +
  labs(title = "Posterior Probability per Sample (allele)",
       x = "Posterior Probability",
       y = "Density") +
  facet_wrap(~patient_id, scales = "free_y") +
  theme_minimal()
ggsave("8.6_allele_pcnv.png", width = 12, height = 8)

# p_cnv
ggplot(all_pcnv, aes(x = p_cnv, fill = patient_id)) +
  geom_density(alpha = 0.5) +
  labs(title = "Posterior Probability per Sample (p_cnv)",
       x = "Posterior Probability",
       y = "Density") +
  facet_wrap(~patient_id, scales = "free_y") +
  theme_minimal()

ggsave("8.6_posterior_pcnv.png", width = 12, height = 8)

# Merge with Seurat data
# Ensure barcodes are standardized
 all_pcnv$barcode <- paste0(all_pcnv$patient_id, "_", all_pcnv$cell)

 # Define function to make the barcodes unique
 make_unique <- function(x) {
   make.unique(x, sep = "__")
 }
# Ensure barcodes are standardized in seurat metadata
barcode_df <- as.data.frame(seu_combined@meta.data)
barcode_df$old_barcode <- rownames(barcode_df)
barcode_df$patient_id <- barcode_df$patient_id %||% barcode_df$Sample %||% NA
new_barcodes <- paste0(barcode_df$patient_id, "_", sub(".*_", "", barcode_df$old_barcode))
new_barcodes <- make_unique(new_barcodes)
colnames(seu_combined) <- new_barcodes
rownames(seu_combined@meta.data) <- new_barcodes
seu_combined$barcode <- rownames(seu_combined@meta.data)

# Filter Seurat to only cells with matching barcodes
matching_barcodes <- intersect(rownames(seu_combined@meta.data), all_pcnv$barcode)
seu_combined <- subset(seu_combined, cells = matching_barcodes)

# Merge metadata
seu_combined@meta.data <- left_join(seu_combined@meta.data, all_pcnv, by = c("barcode" = "barcode"))




#Plot the UMAP with probability
p1 <- ggplot(seu_combined@meta.data, aes(x = UMAP_1, y = UMAP_2, color = p_cnv_y)) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_color_gradient(low = "lightgrey", high = "firebrick") +
  labs(title = "Posterior CNV Probability (p_cnv_y)", color = "CNV\nProbability") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    strip.text = element_text(size = 14),
    panel.grid = element_blank(),
    strip.background = element_blank(), 
    legend.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5))

pdf("8.6_UMAP_cnvprob.pdf", width = 7, height = 7)  
p1
dev.off()

# Set a confidence treshold
seu_combined@meta.data$cnv_confidence <- ifelse(seu_combined@meta.data$p_cnv_y > 0.9, "High", "Low")

p2 <- ggplot(seu_combined@meta.data, aes(x = UMAP_1, y = UMAP_2, color = cnv_confidence)) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_color_manual(values = c("Low" = "grey", "High" = "darkred")) +
  labs(title = "UMAP by CNV Call Confidence", color = "Confidence>0.9") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    strip.text = element_text(size = 14),
    panel.grid = element_blank(),
    strip.background = element_blank(), 
    legend.position = "right",
    plot.title = element_text(hjust = 0.5))
p2

pdf("8.6_UMAP_cnvconfidence1.pdf", width = 7, height = 7)  
p2
dev.off()

# Same plot with only tumor cells (it would be good to show that T cell that were identified as tumor cells has low confidence level)

p3 <- seu_combined@meta.data  %>% 
  filter(compartment_opt.x=="tumor") %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cnv_confidence)) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_color_manual(values = c("Low" = "grey", "High" = "darkred")) +
  labs(title = "Tumor cells  by CNV Call Confidence", color = "Confidence>0.9") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    strip.text = element_text(size = 14),
    panel.grid = element_blank(),
    strip.background = element_blank(), 
    legend.position = "right",
    plot.title = element_text(hjust = 0.5))
p3

pdf("8.6_UMAP_tumor_cells_cnvconfidence.pdf", width = 7, height = 7)  
p3
dev.off()
