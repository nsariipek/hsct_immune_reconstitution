# Nurefsan Sariipek, 250220- Updated at 250428
# Visualize the souporcell results 
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

#Load the saved souporcell result table
final_dataset <- read_csv("~/250428_final_dataset.csv") 

# For Peter:
#final_dataset <- read_csv("AuxiliaryFiles/final_dataset_250418.csv")

#reorder
celltypes <-  c("HSC MPP", "MEP", "EoBasoMast Precursor", "Megakaryocyte Precursor", "Early Erythroid", "Late Erythroid", 
                "LMPP", "Cycling Progenitor", "Early GMP", "Late GMP", "Pro-Monocyte", "CD14 Mono", 
                "CD16 Mono", "cDC", "pDC", "Early Lymphoid", "Pro-B", "Pre-B", 
                "Immature B", "Mature B", "Plasma Cell", "CD4 Naive", "CD4 Central Memory", "CD4 Effector Memory", 
                "CD4 Regulatory", "CD8 Naive", "CD8 Central Memory", "CD8 Effector Memory 1", "CD8 Effector Memory 2", 
                "CD8 Tissue Resident Memory", "T Proliferating", "NK", "NK CD56high", "NK Proliferating", "Stromal")

# Convert to factor with specified order
final_dataset$celltype <- factor(final_dataset$celltype, levels = celltypes)

# Define the souporcell colors
souporcell_colors <-  c("donor" = "#4B3140",recipient ="#E4C9B0", "unknown" = "#b0b0b0")
# cohort colors
cohort_colors <- c("long-term-remission" = "#546fb5FF","relapse" = "#e54c35ff")

# Organize the dataset 
t1 <- final_dataset %>%
  drop_na(celltype) %>%
  filter(sample_status=="remission") %>%
  count(patient_id, celltype, sample_status, origin, name = "count") %>%
  group_by(patient_id, celltype) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

p1 <- t1 %>%
  filter(origin %in% c("donor", "recipient")) %>%
  ggplot(aes(x = patient_id, y = proportion, fill = origin)) +
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
    panel.grid = element_blank(),                     # Remove grids
    axis.line = element_line(color = "black", size = 0.5),  # Keep normal axis lines
    axis.ticks = element_line(color = "black", size = 0.5), # Only ticks thinner
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
pdf("6.4_souporcell_results_all_cells_.pdf", width = 16, height = 8)
p1
dev.off()
# See the genotype results per patient by each sample status for the ones that worked
status <- c("pre-transplant","remission","relapse")
final_dataset$sample_status <- factor(final_dataset$sample_status, levels = status)

t2 <- final_dataset %>%
  drop_na(celltype) %>%
  filter(origin %in% c("donor", "recipient")) %>%
  count(patient_id,origin,sample_status, name = "count") %>%
  group_by(patient_id, sample_status) %>%  
  mutate(proportion = count / sum(count)) %>%
  ungroup()

p2 <- ggplot(t2, aes(x = sample_status, y = proportion, fill = origin)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ patient_id) +
  scale_fill_manual(values = souporcell_colors) +
  labs(x = "Patient ID", y = "Proportion", fill = "Cell Origin",
       title = "Proportion of Cell Origins (Donor vs. Recipient) in Each Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2

# Save as a pdf file 
pdf("6.4_souporcell_per_patient.pdf", width = 12, height = 8)
p2
dev.off()

#### Heatmap ####
# Prepare data
t3 <- final_dataset %>%
  drop_na(celltype) %>%
 filter(origin %in% c("donor", "recipient"),
    sample_status == "remission",
    timepoint %in% c("3","5","6")) %>%
  count(patient_id, celltype, origin, cohort, TP53_status, name = "count") %>%
  group_by(patient_id, celltype, cohort, TP53_status) %>%
  mutate(
    total_cells = sum(count),
    donor_percentage = sum(ifelse(origin == "donor", count, 0)) / total_cells * 100
  ) %>%
  ungroup() %>%
  filter(total_cells > 0) %>%
  select(patient_id, cohort, celltype, donor_percentage, TP53_status) %>%
  distinct()

# Save the table to put intocelltypentary
write_csv(t3,"donorchimerism.csv")


# Plot
plot <- t3 %>% 
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
plot

pdf("6.4_souporcell_heatmap.pdf", width = 8, height = 8)
heatmap
dev.off()

# Compare the donor chimerism ratio between cohorts in 100-180 days samples per each celltype
donor_chimerism_comparision <- t3 %>%
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

donor_chimerism_comparision
pdf("6.4_donor_chimerism_comparision_all.pdf", width = 8, height = 6)
donor_chimerism_comparision
dev.off()

# Just progenitor cells
plot_prog <- function(df, title = NULL) {
  df %>%
    filter(celltype %in% c( "HSC MPP","MEP","LMPP","Cycling Progenitors","Early GMP","Late GMP")) %>%
    ggplot(aes(x = cohort, y = donor_percentage, color = cohort)) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) +
    facet_wrap(~ celltype, nrow = 1) +
    stat_compare_means(aes(group = cohort), method = "wilcox.test",
                       label.y = 50, label.x = 1.25, size = 3, label = "p.format") +
    labs(y = "Donor Chimerism", title = title) +
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
}

# Generate both plots
prog_all <- plot_prog(t3, title = "All Cells")
prog_MT  <- plot_prog(filter(t3, TP53_status == "MUT"), title = "TP53-Mutated Cells")

pdf("6.4_donor_chimerism_progenitors_all.pdf", width = 8, height = 6)
prog_all 
dev.off()

pdf("6.4_donor_chimerism_progenitors_TP53.pdf", width = 8, height = 6)
prog_MT 
dev.off()



