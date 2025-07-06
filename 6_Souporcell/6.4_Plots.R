# Nurefsan Sariipek and Peter van Galen, 250702
# Generate dot plots of donor/recipient proportions

# Load the libraries
library(tidyverse)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(janitor)
library(ggtext)

# Empty environment
rm(list=ls())

# For Nurefsan:
setwd("~/TP53_ImmuneEscape/6_Souporcell/")
# For Peter:
#setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/6_Souporcell/")

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Define the souporcell colors
souporcell_colors <-  c("donor" = "#4B3140", recipient ="#E4C9B0", unknown = "#b0b0b0")

# Cohort colors
cohort_colors <- c("long-term-remission" = "#546fb5", relapse = "#e54c35")

# Calculate the proportion recipient and donor cells
t1 <- as_tibble(seu@meta.data) %>%
  filter(!is.na(celltype), sample_status == "remission",
         souporcell_origin %in% c("donor", "recipient")) %>%
  mutate(souporcell_origin = factor(souporcell_origin, levels = c("donor", "recipient"))) %>%
  group_by(patient_id, celltype, sample_status, souporcell_origin) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient_id, celltype) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

p1 <- t1 %>%
  ggplot(aes(x = patient_id, y = proportion, fill = souporcell_origin)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  facet_wrap(~celltype, scales = "free_x", nrow = 3) +
  scale_fill_manual(values = souporcell_colors) +
  scale_x_discrete(labels = function(x) substr(x, 1, 3), expand = c(0.05, 0.05)) +
  labs(
    x = "Patient ID",
    y = "Proportion",
    fill = "Genotype") +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title.x = element_text(size = 6, color = "black"),
    axis.title.y = element_text(size = 6, color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    strip.text = element_text(size = 8, color = "black"),
    panel.spacing = unit(1.5, "lines"),
    plot.title = element_text(size = 11, color = "black", hjust = 0.5),
    legend.text = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 8, color = "black"),
    legend.position = "right")

# View the plot
p1

# Save as a pdf file 
pdf("6.4.1_Souporcell_all_barplots.pdf", width = 16, height = 8)
p1
dev.off()

# See the genotype results per patient by each sample status for the ones that worked
t2 <- as_tibble(seu@meta.data) %>%
  filter(!is.na(celltype),
         souporcell_origin %in% c("donor", "recipient")) %>%
  mutate(souporcell_origin = factor(souporcell_origin, levels = c("donor", "recipient"))) %>%
  group_by(patient_id, souporcell_origin, sample_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient_id, sample_status) %>%  
  mutate(proportion = count / sum(count)) %>%
  ungroup()

p2 <- ggplot(t2, aes(x = sample_status, y = proportion, fill = souporcell_origin)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ patient_id) +
  scale_fill_manual(values = souporcell_colors) +
  labs(x = "Patient ID", y = "Proportion", fill = "Cell Origin",
       title = "Proportion of cell origins in each sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2

# Save as a pdf file 
pdf("6.4.2_Souporcell_per_patient.pdf", width = 12, height = 8)
p2
dev.off()

# Prepare data for heatmap visualization
t3 <- as_tibble(seu@meta.data) %>%
 filter(!is.na(celltype),
        souporcell_origin %in% c("donor", "recipient"),
        sample_status == "remission",
        timepoint %in% c("3","5","6")) %>%
  group_by(patient_id, celltype, souporcell_origin, cohort, TP53_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient_id, celltype, cohort, TP53_status) %>%
  mutate(
    total_cells = sum(count),
    donor_percentage = sum(ifelse(souporcell_origin == "donor", count, 0)) / total_cells * 100
  ) %>%
  ungroup() %>%
  filter(total_cells > 0) %>%
  select(patient_id, cohort, celltype, donor_percentage, TP53_status) %>%
  distinct()

# Plot
heatmap <- t3 %>% 
ggplot(aes(x = celltype, y = patient_id, fill = donor_percentage)) +
 geom_raster()+
  scale_fill_gradientn(
    colors = c("#E4C9B0", "#C9AB8F", "#AA8D6E", "#866A78", "#5F4B5B", "#4B3140"),
    limits = c(0, 100),
    name = "Percentage",
    na.value = "grey50"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev(unique(t3$patient_id))) +
  labs(
    x = "Cell Type",
    y = "Patient ID",
    title = "Donor Chimerism by Souporcell in 3-6 mo Remission Samples") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    plot.title = element_text(size = 12, color = "black", hjust = 0.5),  
    legend.title = element_text(size = 11, color = "black"),
    panel.grid = element_blank(),                
    axis.ticks = element_line(color = "black", size = 0.5)) +
    coord_fixed(ratio = 1)  # Make tiles square

# View the heatmap
heatmap

pdf("6.4.3_Souporcell_heatmap.pdf", width = 8, height = 8)
heatmap
dev.off()

# Visualize the donor chimerism ratio between cohorts 3-6M samples for each celltype
donor_chimerism_comparison <- t3 %>%
  ggplot(aes(x = cohort, y = donor_percentage, color = cohort)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +  
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) + 
  facet_wrap(~ celltype, nrow = 5) +
  stat_compare_means(aes(group = cohort), method = "wilcox.test", label.y = 50, label.x = 1.25, size= 3, label="p.format") +  # Mann-Whitney U
  theme_bw() +
  labs( x = "Cohort",
       y = "Donor Chimerism") +
  scale_color_manual(values = cohort_colors) +
  scale_y_continuous(limits = c(0, 100)) + 
  theme_minimal(base_size = 10) +
  theme(#aspect.ratio = 2,
    axis.line = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 10, color = "black"),
   legend.title = element_text(size = 11, color = "black"),
    panel.grid = element_blank(),                
    axis.ticks = element_line(color = "black", size = 0.5)) 

# View the plot
donor_chimerism_comparison

pdf("6.4.4_Jitter_chimerism_all.pdf", width = 8, height = 6)
donor_chimerism_comparison
dev.off()


# Now for progenitors only
t4 <- as_tibble(seu@meta.data) %>%
  filter(!is.na(celltype),
         souporcell_origin %in% c("donor", "recipient"),
         sample_status == "remission",
         timepoint %in% c("3","5","6")) %>%
  filter(celltype %in% c("HSC MPP","MEP","LMPP","Cycling Progenitors","Early GMP"))

t4$celltype_merged <- "Merged_Progenitors"

merged_per_tb <- t4 %>%
group_by(patient_id, celltype_merged, souporcell_origin, cohort, TP53_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient_id, celltype_merged, cohort, TP53_status) %>%
  mutate(
    total_cells = sum(count),
    donor_percentage = sum(ifelse(souporcell_origin == "donor", count, 0)) / total_cells * 100
  ) %>%
  ungroup() %>%
  filter(total_cells > 0) %>%
  select(patient_id, cohort, celltype_merged, donor_percentage, TP53_status) %>%
  distinct()

merged_prog <- merged_per_tb %>% 
  ggplot(aes(x = cohort, y = donor_percentage, color = cohort)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +  
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) +
  stat_compare_means(aes(group = cohort), method = "wilcox.test", label.y = 50, label.x = 1.25, size= 3, label="p.format") +  # Mann-Whitney U
  theme_bw() +
  labs(x = "Cohort",
       y = "Donor Chimerism") +
  scale_color_manual(values = cohort_colors) +
  scale_y_continuous(limits = c(0, 101)) + 
  theme_minimal(base_size = 10) +
  theme( aspect.ratio = 2,
         axis.line = element_line(color = "black", size = 0.5),
         axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, color = "black"),
         axis.text.y = element_text(size = 8, color = "black"),
         axis.title.x = element_text(size = 10, color = "black"),
         axis.title.y = element_text(size = 11, color = "black"),
         strip.text = element_text(size = 10, color = "black"),
         legend.title = element_text(size = 11, color = "black"),
         panel.grid = element_blank(),                
         axis.ticks = element_line(color = "black", size = 0.5)) 

# View the plot
merged_prog

pdf("6.4.5_Donor_chimerism_merged_HSPC.pdf", width = 8, height = 6)
merged_prog
dev.off()

# Just progenitor cells
prog_all <- t3 %>%
    filter(celltype %in% c( "HSC MPP","MEP","LMPP","Cycling Progenitor","Early GMP")) %>%
    ggplot(aes(x = cohort, y = donor_percentage, color = cohort)) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) +
    facet_wrap(~ celltype, nrow = 1) +
    stat_compare_means(aes(group = cohort), method = "wilcox.test",
                       label.y = 50, label.x = 1.25, size = 3, label = "p.format") +
    labs(y = "Donor Chimerism") +
    scale_color_manual(values = cohort_colors) +
    coord_cartesian(ylim = c(0, 100)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.line = element_line(color = "black", size = 0.5),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.x = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 11, color = "black"),
      strip.text = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 11, color = "black"),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.5), 
      aspect.ratio = 1
    )

# View
prog_all

pdf("6.4.6_Donor_chimerism_progenitors_types.pdf", width = 8, height = 6)
prog_all
dev.off()


# All cells
t5 <- as_tibble(seu@meta.data) %>%
  filter(!is.na(celltype),
         souporcell_origin %in% c("donor", "recipient"),
         sample_status == "remission",
         timepoint %in% c("3","5","6")) %>%
group_by(patient_id, souporcell_origin, cohort) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient_id, cohort) %>%
  mutate(
    total_cells = sum(count),
    donor_percentage = sum(ifelse(souporcell_origin == "donor", count, 0)) / total_cells * 100
  ) %>%
  ungroup() %>%
  filter(total_cells > 0) %>%
  select(patient_id, cohort, donor_percentage) %>%
  distinct()

merged_all <- t5 %>% 
  ggplot(aes(x = cohort, y = donor_percentage, color = cohort)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +  
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) +
  stat_compare_means(aes(group = cohort), method = "wilcox.test", label.y = 50, label.x = 1.25, size= 3, label="p.format") +  # Mann-Whitney U
  theme_bw() +
  labs( x = "Cohort",
        y = "Donor Chimerism") +
  scale_color_manual(values = cohort_colors) +
  scale_y_continuous(limits = c(0, 101)) + 
  theme_minimal(base_size = 10) +
  theme( aspect.ratio = 2,
         axis.line = element_line(color = "black", size = 0.5),
         axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, color = "black"),
         axis.text.y = element_text(size = 8, color = "black"),
         axis.title.x = element_text(size = 10, color = "black"),
         axis.title.y = element_text(size = 11, color = "black"),
         strip.text = element_text(size = 10, color = "black"),
         legend.title = element_text(size = 11, color = "black"),
         panel.grid = element_blank(),                
         axis.ticks = element_line(color = "black", size = 0.5)) 

pdf("6.4.7_Donor_chimerism_merged_all.pdf", width = 8, height = 6)
merged_all
dev.off()
