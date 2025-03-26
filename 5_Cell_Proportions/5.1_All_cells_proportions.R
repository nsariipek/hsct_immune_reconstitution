# Calculate the myeloid cells proportions
# Load the libraries
library(tidyverse)
library(janitor)
library(ggrepel)
library(ggpubr)

# Empty environment
rm(list=ls())

# Set directory
setwd("~/TP53_ImmuneEscape/5_Cell_Proportions/")

# # Load the saved seurat object and turn into a df
# seu <- readRDS("~/250128_seurat_annotated_final.rds")
# 
# # turn into a metadata
# seu_df <- seu@meta.data 
# 
# #change the wrong sample identification
# seu_df <- seu_df %>%
#   mutate(sample_status = if_else(patient_id == "P32" & sample_status == "remission", 
#                                  "relapse", sample_status))
# #save the metadata for future use 
# write.csv(seu_df, "~/seu_df_250325.csv", row.names = FALSE)

# Load the seurat df
seu_df <- read_csv("~/seu_df_250325.csv")


celltype_levels <- c(
  "Progenitors", "Early Erythroids", "Mid Erythroids", "Late Erythroids",
  "Pro Monocytes", "Monocytes", "Non Classical Monocytes", "cDC", "pDC",
  "Pro B cells", "Pre-B", "B cells", "Plasma cells",
  "CD4 Naïve", "CD4 Effector Memory", "CD4 Memory", "Treg",
  "CD8 Naïve", "CD8 Effector", "CD8 Memory", "CD8 Exhausted",
  "γδ T", "NK T", "Adaptive NK", "CD56 Bright NK", "CD56 Dim NK",
  "Cycling T-NK cells", "UD1", "UD2", "UD3"
)

# Convert celltype column to a factor with desired order
seu_df$celltype <- factor(seu_df$celltype, levels = celltype_levels)


# Set the sample order same as Figure 1 swimmer plot
sample_order <- c(paste0("P", str_pad(1:4, 2, pad = "0")),     
                  paste0("P", str_pad(13:27, 2, pad = "0")),     
                  paste0("P", str_pad(5:12, 2, pad = "0")),   
                  paste0("P", str_pad(28:33, 2, pad = "0")))     

bar_data <- seu_df %>%
  mutate(patient_id = factor(patient_id, levels = sample_order)) %>%
  filter(library_type == "MNC",timepoint %in% c("3","5","6") & sample_status =="remission" )%>%
  group_by(patient_id, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(percent = count / sum(count) * 100)

# Set the color scale

celltype_colors <- c(
  "Progenitors" = "#3B1B53FF",
  "Early Erythroids" = "#D60047FF",
  "Mid Erythroids" = "#924822FF",
  "Late Erythroids" = "#AE1F63FF",
  "Pro Monocytes" = "#99CC00FF",
  "Monocytes" = "#E4AF69FF",
  "Non Classical Monocytes" = "#7A65A5FF",
  "cDC" = "#5DB1DDFF",
  "pDC" = "#CDDEB7FF",
  "Pro B cells" = "#14FFB1FF",
  "Pre-B" = "#00991AFF",
  "B cells" = "#003399FF",
  "Plasma cells" = "#802268FF",
  "CD4 Naïve" = "#466983FF",
  "CD4 Effector Memory" = "#D58F5CFF",
  "CD4 Memory" = "#C75127FF",
  "Treg" = "#FFC20AFF",
  "CD8 Naïve" = "#33CC00FF",
  "CD8 Effector" = "#612A79FF",
  "CD8 Memory" = "#0099CCFF",
  "CD8 Exhausted" = "#CE3D32FF",
  "γδ T" = "#D595A7FF",
  "NK T" = "#5050FFFF",
  "Adaptive NK" = "#1A0099FF",
  "CD56 Bright NK" = "#00D68FFF",
  "CD56 Dim NK" = "#008099FF",
  "Cycling T-NK cells" = "#F0E685FF",
  "UD1" = "#A9A9A9FF",
  "UD2" = "#837B8DFF",
  "UD3" = "#5A655EFF")

p1 <- ggplot(bar_data, aes(x = patient_id, y = percent, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "Patient ID", y = "Cell Type Proportion (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )
p1

# Save as a pdf
pdf("5.1_Allcells_proportions_3-6mo.pdf", width = 10, height = 5)
p1
dev.off()

################################################################################################################

# Calculate total cells and myeloid cells per patient ==percentage of all myeloid cells out of all cell types in each sample (MNC libraries only)
# 
# myeloid_types <- c("Progenitors","Early Erythroids","Mid Erythroids",
#                    "Late Erythroids","pDC","cDC","Pro Monocytes", 
#                    "Monocytes","Non Classical Monocytes")

proportions_df <- seu_df %>%
  filter(timepoint %in% c("3", "5", "6"),
         sample_status == "remission", library_type=="MNC") %>%
  group_by(patient_id, survival) %>%
  mutate(total_cells = n()) %>%  # total cells per sample
 # filter(celltype %in% myeloid_types) %>%
  count(patient_id, survival, total_cells, celltype, name = "cell_count") %>%
  mutate(percent = (cell_count / total_cells) * 100) %>%
mutate(survival = factor(survival, levels = c("Non-relapsed", "Relapsed")))
# %>% 
#   mutate(celltype = factor(celltype, levels = myeloid_types))

p2 <- ggplot(proportions_df, aes(x = survival, y = percent, fill = survival)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(shape = 21, size = 2, stroke = 0.2, color = "black", width = 0.15) +
  facet_wrap(~ celltype, ncol = 10) +
  scale_fill_manual(values = c("Relapsed" = "#E64B35FF", "Non-relapsed" = "#4DBBD5FF")) +
  labs(y = "% of all cells", x = NULL) +
  ggpubr::stat_compare_means(
    aes(x = survival, y = percent, group = survival),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    label.y = 40
  ) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 15, face = "plain", color = "black"),
    strip.text = element_text(size = 10, face = "plain", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "none",
    aspect.ratio = 2,
    plot.margin = margin(4, 4, 4, 4)
  )

p2
# Save as a pdf
pdf("5.1_Allcell_per_celltype_proportions_3-6mo.pdf", width = 10, height = 8)
p2
dev.off()

# Addition if you want to calculate the p-values separately 

library(broom)

pvals_df <- proportions_df %>%
  group_by(celltype) %>%
  filter(!is.na(percent), !is.na(survival)) %>%
  summarise(
    p.value = tryCatch(wilcox.test(percent ~ survival)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(label = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

