# Nurefsan Sariipek, 231103 updated at 250225
# Check T/NK cell proportions
# Load the libraries
library(tidyverse)
library(janitor)
library(RColorBrewer)
library(randomcoloR)
library(ggrepel)
library(ggpubr)
library(ggnewscale)

# Empty environment
rm(list=ls())

# For Nurefsan:
setwd("~/TP53_ImmuneEscape/5_Cell_Proportions/")

# Load the seurat df
seu_df <- read_csv("~/seu_df_250325.csv")
seu_df <- seu_df %>%
  mutate(TP53_status = case_when(
    patient_id %in% c("P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P12", "P14", "P17") ~ "MT",
    patient_id %in% c("P13", "P15", "P16", "P18", "P19", "P20", "P21", "P22", "P23", "P24", "P25", "P26", "P27",
                      "P28", "P29", "P30", "P31", "P32", "P33") ~ "WT",
    TRUE ~ NA_character_  # fallback in case of unmatched ID
  ))


# Set the sample order same as Figure 1 swimmer plot
sample_order <- c(paste0("P", str_pad(1:4, 2, pad = "0")),     
                  paste0("P", str_pad(13:27, 2, pad = "0")),     
                  paste0("P", str_pad(5:12, 2, pad = "0")),   
                  paste0("P", str_pad(28:33, 2, pad = "0")))     

cell_order <- c("CD4 Naïve", "CD4 Memory", "CD4 Effector Memory", "CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Exhausted",
  "Treg","γδ T")


bar_data <- seu_df %>%
  mutate(patient_id = factor(patient_id, levels = sample_order)) %>%
  filter(timepoint %in% c("3","5","6") & sample_status =="remission")%>%
  filter(celltype %in% c("CD4 Naïve","CD4 Memory","CD4 Effector Memory","Treg","CD8 Naïve","CD8 Memory",
                         "CD8 Effector","CD8 Exhausted", "γδ T")) %>%
  group_by(patient_id, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(percent = count / sum(count) * 100) %>%
  

 bar_data$celltype <- factor(bar_data$celltype, levels = cell_order)


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
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  )
p1
# Save as a pdf
pdf("5.2_Tcells_proportions_3-6mo.pdf", width = 10, height = 8)
p1
dev.off()

###############################################################
# Calculate the proportions for T cells

  proportions_df <- seu_df %>%
  filter(celltype %in% c("CD4 Naïve","CD4 Memory","CD4 Effector Memory","Treg",
                         "CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Exhausted", 
                         "delta-gamma T") &
           timepoint %in% c("3","5","6") &
           sample_status == "remission"
        # &  TP53_status == "MT"
         ) %>%
  group_by(sample_id, survival) %>%
  reframe(tabyl(celltype)) %>%
  group_by(sample_id, survival) %>%
  mutate(
    total_T_cells = sum(n),
    percent_within_T = (n / total_T_cells) * 100,
    celltype = factor(celltype, levels = c("CD4 Naïve","CD4 Memory","CD4 Effector Memory",
                                           "Treg","CD8 Naïve","CD8 Memory",
                                           "CD8 Effector","CD8 Exhausted","delta-gamma T")),
    survival = factor(survival, levels = c("Non-relapsed", "Relapsed"))
  )



# Plot
p2 <- ggplot(proportions_df, aes(x = survival, y = percent_within_T, fill = survival)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9, linewidth = 0.5, color = "black") + 
  geom_jitter(shape = 21, size = 1.8, width = 0.15, stroke = 0.2, color = "black", aes(fill = survival)) +
  coord_cartesian(ylim = c(0, 50)) +
  facet_wrap(~ celltype, ncol = 10) +
  scale_fill_manual(values = c("Relapsed" = "#E64B35FF", "Non-relapsed" = "#4DBBD5FF")) +
  labs(y = "% within total T cells", x = NULL) +
  stat_compare_means(
    aes(x = survival, y = percent_within_T, group = survival),  # safest
    method = "wilcox.test",
    label = "p.format",
    size = 2.3,
    hide.ns = TRUE
  ) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title.y = element_text(size = 12, face = "plain"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "none",
    aspect.ratio = 2,
    plot.margin = margin(4, 4, 4, 4)
  )


# Check the plot
p2
# Save as a pdf
pdf("5.2_T_proportions_TP53_MT_value.pdf", width = 8, height = 6)
p2
dev.off()




