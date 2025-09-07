# Nurefsan Sariipek and Peter van Galen, 250706
# Generate plots of donor/recipient proportions

# Load the libraries
library(tidyverse)
library(Seurat)
library(ggpubr)
library(ggtext)

# Empty environment
rm(list = ls())

# For VM:
# setwd("~/hsct_immune_reconstitution/08_Souporcell/")
# For Peter:
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/08_Souporcell/")

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# SET UP -----------------------------------------------------------------

# Extract & filter metadata
metadata_tib <- as_tibble(seu@meta.data) %>%
  filter(
    !is.na(celltype),
    souporcell_origin %in% c("donor", "recipient")
  ) %>%
  mutate(
    souporcell_origin = factor(
      souporcell_origin,
      levels = c("donor", "recipient")
    )
  )

# Define the souporcell colors
souporcell_colors <- c(
  donor = "#4B3140",
  recipient = "#E4C9B0",
  unknown = "#b0b0b0"
)

# Cohort colors
cohort_colors <- c("long-term-remission" = "#546fb5", relapse = "#e54c35")


# 1. BARPLOTS -----------------------------------------------------------------

# Chimerism per patient by sample status
t1 <- metadata_tib %>%
  group_by(patient_id, souporcell_origin, sample_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient_id, sample_status) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Barplot per patient by sample status
p1 <- t1 %>%
  ggplot(aes(x = sample_status, y = proportion, fill = souporcell_origin)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~patient_id, ncol = 7) +
  scale_fill_manual(values = souporcell_colors) +
  labs(
    x = "Patient ID",
    y = "Proportion",
    fill = "Cell origin",
    title = "Proportion of cell origins in each sample"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# View
p1

# Save as a pdf file
pdf("8.4.1_Barplots_per_patient.pdf", width = 8, height = 8)
p1
dev.off()

# Chimerism per cell type by patient
t2 <- metadata_tib %>%
  filter(sample_status == "remission", timepoint %in% c(3, 5, 6)) %>%
  select(patient_id, celltype, souporcell_origin) %>%
  # Add cell type percent
  arrange(celltype) %>%
  add_count(celltype, name = "celltype_count") %>%
  mutate(celltype_percent = celltype_count / n() * 100) %>%
  mutate(
    celltype = paste0(celltype, " (", round(celltype_percent, 2), "%)")
  ) %>%
  mutate(
    celltype = factor(
      celltype,
      levels = unique(celltype)
    )
  ) %>%
  # Add souporcell donor/recipient proportion
  group_by(
    patient_id,
    celltype,
    souporcell_origin
  ) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient_id, celltype) %>%
  mutate(souporcell_proportion = count / sum(count))

# Barplot per cell type by patient
p2 <- t2 %>%
  ggplot(aes(
    x = patient_id,
    y = souporcell_proportion,
    fill = souporcell_origin
  )) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  facet_wrap(~celltype, scales = "free_x", nrow = 5) +
  scale_fill_manual(values = souporcell_colors) +
  scale_x_discrete(
    labels = function(x) substr(x, 1, 3),
    expand = c(0.05, 0.05)
  ) +
  labs(
    x = "Patient ID",
    y = "Proportion",
    fill = "Genotype"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, , vjust = 0.5, hjust = 1),
    axis.ticks = element_line(color = "black"),
    panel.spacing = unit(1.5, "lines"),
    plot.title = element_text(size = 10, color = "black", hjust = 0.5),
  )

# View
p2

# Save as a pdf file
pdf("8.4.2_Barplots_per_celltype.pdf", width = 16, height = 12)
p2
dev.off()


# 2. HEATMAP AND JITTER PLOTS -------------------------------------------------

# Prepare data for heatmap visualization
counts_tib <- metadata_tib %>%
  filter(sample_status == "remission", timepoint %in% c(3, 5, 6)) %>%
  dplyr::count(patient_id, cohort, celltype, souporcell_origin)

# Pivot wider to compute proportion donor
props_tib <- counts_tib %>%
  pivot_wider(
    names_from = souporcell_origin,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(donor_percentage = donor / (donor + recipient) * 100)

# For paper text:
props_tib %>% filter(celltype == "Stromal")
props_tib %>% filter(celltype == "Plasma Cell") %>% pull(donor) %>% sum
props_tib %>% filter(celltype == "Plasma Cell") %>% pull(recipient) %>% sum
649 / (649 + 451)
# "The replacement of recipient cell by donor cells varied across cell populations, with the persistence of recipient cells being highest for stromal cells and plasma cells (100% and 59%, respectively)"

# HEATMAP ---------------------------------------

heatmap <- props_tib %>%
  mutate(patient_id = factor(patient_id, levels = unique(patient_id))) %>%
  complete(patient_id, celltype) %>%
  ggplot(aes(x = celltype, y = patient_id, fill = donor_percentage)) +
  geom_raster() +
  scale_fill_gradientn(
    colors = c("#E4C9B0", "#4B3140"),
    limits = c(0, 100),
    name = "Donor percentage",
    na.value = "grey80"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(
    expand = c(0, 0),
    limits = rev(unique(props_tib$patient_id))
  ) +
  labs(
    y = "Patient ID",
    title = "Donor chimerism in 3-6 mo remission samples"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    axis.ticks = element_line(color = "black")
  ) +
  coord_fixed(ratio = 1) # Make tiles square

# View
heatmap

pdf("8.4.3_Souporcell_heatmap.pdf", width = 9, height = 5)
heatmap
dev.off()

# JITTER ALL ------------------------------------

# Donor chimerism per celltype by cohort
jitter_all_plot <- props_tib %>%
  ggplot(aes(x = cohort, y = donor_percentage, color = cohort)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    size = 0.2,
    color = "#00000080"
  ) +
  facet_wrap(~celltype, nrow = 5) +
  stat_compare_means(
    aes(group = cohort),
    method = "wilcox.test",
    label.y = 50,
    label.x = 0.8,
    size = 3,
    label = "p.format",
    show.legend = F
  ) +
  labs(x = "Cohort", y = "Donor chimerism") +
  scale_color_manual(values = cohort_colors) +
  scale_y_continuous(limits = c(-1, 102)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank()
  )

# View the plot. Warning is due to Stromal stats
jitter_all_plot

pdf("8.4.4_Jitter_all.pdf", width = 8, height = 6)
jitter_all_plot
dev.off()

# JITTER MERGED ALL -----------------------------

# Donor chimerism per patient (cell types merged)
merged_props_tib <- props_tib %>%
  group_by(patient_id, cohort) %>%
  summarize(donor = sum(donor), recipient = sum(recipient)) %>%
  mutate(donor_percentage = donor / (donor + recipient) * 100)

# Plot
merged_all_plot <- merged_props_tib %>%
  ggplot(aes(x = cohort, y = donor_percentage, color = cohort)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
  coord_cartesian(ylim = c(0, 100)) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    size = 0.4,
    color = "#00000080"
  ) +
  stat_compare_means(
    aes(group = cohort),
    method = "wilcox.test",
    label.y = 50,
    label.x = 1.25,
    size = 3,
    label = "p.format",
    show.legend = F
  ) +
  labs(y = "Donor chimerism") +
  scale_color_manual(values = cohort_colors) +
  theme_bw() +
  theme(
    aspect.ratio = 2,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank()
  )

# View
merged_all_plot

pdf("8.4.5_Jitter_all_merged.pdf", width = 4, height = 4)
merged_all_plot
dev.off()


# JITTER HSPCs ----------------------------------

# Jitter plot per HSPC cell type by cohort
prog_plots <- props_tib %>%
  filter(
    celltype %in% c("HSC MPP", "MEP", "LMPP", "Cycling Progenitor", "Early GMP")
  ) %>%
  ggplot(aes(x = cohort, y = donor_percentage, color = cohort)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    size = 0.2,
    color = "#00000080",
    show.legend = F
  ) +
  facet_wrap(~celltype, nrow = 1) +
  stat_compare_means(
    aes(group = cohort),
    method = "wilcox.test",
    label.y = 50,
    label.x = 1.25,
    size = 3,
    label = "p.format",
    show.legend = F
  ) +
  labs(y = "Donor chimerism") +
  scale_color_manual(values = cohort_colors) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  theme(
    aspect.ratio = 2,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank()
  )

# View
prog_plots

pdf("8.4.6_Jitter_HSPCS.pdf", width = 8, height = 4)
prog_plots
dev.off()


# JITTER MERGED HSPCs ---------------------------

# Donor chimerism per patient (cell types merged)
merged_prog_props_tib <- props_tib %>%
  filter(
    celltype %in%
      c("HSC MPP", "MEP", "LMPP", "Cycling Progenitor", "Early GMP")
  ) %>%
  group_by(patient_id, cohort) %>%
  summarize(donor = sum(donor), recipient = sum(recipient)) %>%
  mutate(donor_percentage = donor / (donor + recipient) * 100)

merged_prog_plot <- merged_prog_props_tib %>%
  ggplot(aes(x = cohort, y = donor_percentage, color = cohort)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.5) +
  coord_cartesian(ylim = c(0, 101)) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    size = 0.4,
    color = "#00000080"
  ) +
  stat_compare_means(
    aes(group = cohort),
    method = "wilcox.test",
    label.y = 50,
    label.x = 1.25,
    size = 3,
    label = "p.format",
    show.legend = F
  ) +
  labs(y = "Donor chimerism (HSPCs)") +
  scale_color_manual(values = cohort_colors) +
  theme_bw() +
  theme(
    aspect.ratio = 2,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank()
  )

# View
merged_prog_plot

pdf("8.4.7_Jitter_merged_HSPC.pdf", width = 4, height = 4)
merged_prog_plot
dev.off()

# For text:
merged_prog_props_tib %>%
  group_by(cohort) %>%
  summarize(median_chimerism = median(donor_percentage))
# "We found that HSPC donor chimerism was significantly lower in patients who eventually relapsed compared to those who stayed in remission (median 100% and 63%, p=0.0039)"
