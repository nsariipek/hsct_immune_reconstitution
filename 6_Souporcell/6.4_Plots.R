# Nurefsan Sariipek, 250220
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

#Load the saved souporcell result table( This version is corrected for P32)
final_dataset <- read_csv("~/final_dataset_250414.csv")

# For Peter:
#final_dataset <- read_csv("AuxiliaryFiles/final_dataset_250414.csv"
# Reorder for visulization
final_dataset$sample_status <- factor(final_dataset$sample_status, levels = c("pre_transplant","remission","relapse"))

celltypes <- c("Progenitors", "Early Erythroids", "Mid Erythroids", "Late Erythroids", "Pro Monocytes",
               "Monocytes", "Non Classical Monocytes", "cDC",  "pDC", "Pro B cells", "Pre-B", "B cells",
               "Plasma cells", "CD4 Naïve", "CD4 Effector Memory", "CD4 Memory", "Treg", "CD8 Naïve",
               "CD8 Effector", "CD8 Memory", "CD8 Exhausted",  "γδ T", "NK T", "Adaptive NK", "CD56 Bright NK",
               "CD56 Dim NK", "Cycling T-NK cells", "UD1", "UD2", "UD3")

# Convert to factor with specified order
final_dataset$celltype <- factor(final_dataset$celltype, levels = celltypes)

# Set survival order globally
final_dataset$survival <- factor(final_dataset$survival, levels = c("Non-relapsed", "Relapsed"))

# Set global sample_id ordering by survival
ordered_samples <- final_dataset %>%
  distinct(sample_id, survival) %>%
  arrange(survival, sample_id) %>%
  pull(sample_id)

# Apply globally to final_dataset
final_dataset$sample_id <- factor(final_dataset$sample_id, levels = ordered_samples)

# Optional: shortened version for axis labels
final_dataset$sample_id_short <- substr(as.character(final_dataset$sample_id), 1, 3)

# Define the souporcell colors
souporcell_colors <-  c("donor" = "#4B3140",recipient ="#E4C9B0", "unknown" = "#b0b0b0")
# Survival colors
survival_colors <- c("Non-relapsed" = "#546fb5FF","Relapsed" = "#e54c35ff")

# Organize the dataset 
t1 <- final_dataset %>%
  filter(sample_status=="remission") %>%
  filter(celltype %in% c("Mid Erythroids", "Late Erythroids")) %>%
  count(sample_id, celltype, sample_status, origin, name = "count") %>%
  group_by(sample_id, celltype) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Create stacked bar plot-- Erythroid cells
p1 <-  t1 %>%
  ggplot(aes(x = sample_id, y = proportion, fill = origin)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~celltype, scales = "free_x") + # Separate panels for each celltype
  theme_minimal() +
  labs(title = "Genotype Proportions in Erythroid cells in Remission Samples",
       x = "Sample ID",
       y = "Proportion",
       fill = "Genotype") +
  #scale_fill_manual(values = c("0" = "blue", "1" = "red")) + 
  scale_fill_manual(values = souporcell_colors) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8, color= "black"), 
        axis.text.y = element_text(size = 8, color= "black"),
        axis.ticks.x = element_line(), 
        strip.text = element_text(size = 10, color = "black"),   
        plot.title = element_text(size = 11, color = "black", hjust = 0.5),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        panel.spacing = unit(1.5, "lines"), 
        panel.grid = element_blank(),
        legend.position = "right") 


p1
# Save as a pdf file 
pdf("6.4_souporcell_results_erythroids.pdf", width = 12, height = 8)
p1
dev.off()

# Organize the dataset 
t2 <- final_dataset %>%
  filter(sample_status=="remission") %>%
  count(sample_id, celltype, sample_status, origin, name = "count") %>%
  group_by(sample_id, celltype) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()


p2 <- t2 %>%
  filter(!celltype %in% c("UD1", "UD2", "UD3")) %>%
  ggplot(aes(x = sample_id, y = proportion, fill = origin)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  facet_wrap(~celltype, scales = "free_x", nrow = 3) +
  scale_fill_manual(values = souporcell_colors) +
  scale_x_discrete(labels = function(x) substr(x, 1, 3), expand = c(0.05, 0.05)) +
  labs(
    title = "Genotype proportions in 3-6 months remission samples",
    x = "Sample ID",
    y = "Proportion",
    fill = "Genotype"
  ) +
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
    legend.position = "right"
  )


p2
# Save as a pdf file 
pdf("6.4_souporcell_results_all_cells_.pdf", width = 16, height = 8)
p2
dev.off()

t3 <- final_dataset %>%
  filter(origin %in% c("donor", "recipient")) %>%
  count(patient_id,origin,sample_status, name = "count") %>%
  group_by(patient_id, sample_status) %>%  
  mutate(proportion = count / sum(count)) %>%
  ungroup()


p3 <- ggplot(t3, aes(x = sample_status, y = proportion, fill = origin)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ patient_id) +
  scale_fill_manual(values = souporcell_colors) +
  labs(x = "Patient ID", y = "Proportion", fill = "Cell Origin",
       title = "Proportion of Cell Origins (Donor vs. Recipient) in Each Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3

# Save as a pdf file 
pdf("6.4_souporcell_per_patient_.pdf", width = 12, height = 8)
p3
dev.off()

#### Heatmap ####
# Prepare data
t4 <- final_dataset %>%
  filter(origin %in% c("donor", "recipient"),
    sample_status == "remission",
    timepoint %in% c("3","5","6"),
  TP53_status=="MT",
    !celltype %in% c("UD1", "UD2", "UD3")) %>%
  count(sample_id, celltype, origin, survival, name = "count") %>%
  group_by(sample_id, celltype, survival) %>%
  mutate(
    total_cells = sum(count),
    donor_percentage = sum(ifelse(origin == "donor", count, 0)) / total_cells * 100
  ) %>%
  ungroup() %>%
  select(sample_id, celltype, survival, donor_percentage) %>%
  distinct()

write_csv(t4,"donorchimerism.csv")

# Set survival order
t4$survival <- factor(t4$survival, levels = c("Relapsed", "Non-relapsed"))

# Plot
heatmap <- ggplot(t4, aes(x = celltype, y = sample_id, fill = donor_percentage)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#E4C9B0", "#C9AB8F", "#AA8D6E", "#866A78", "#5F4B5B", "#4B3140"),
    limits = c(0, 100),
    name = "Percentage"
  ) +
  scale_y_discrete(limits = rev) +  # Reverse the y-axis
  labs(
    x = "Cell Type",
    y = "Sample ID",
    title = "Donor Chimerism by Souporcell in 3-6 mo Remission Samples"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    plot.title = element_text(size = 12, color = "black", hjust = 0.5),  
    legend.title = element_text(size = 11, color = "black"),
    panel.grid = element_blank(),                
    axis.ticks = element_line(color = "black", size = 0.5) 
    ) +
  coord_fixed(ratio = 1)  # Make tiles square


# Show plot
heatmap

pdf("6.4_souporcell_heatmap.pdf", width = 8, height = 8)
heatmap
dev.off()

# Compare the donor chimerism ratio between cohorts in 100-180 days samples per each celltype
donor_chimerism_comparision <- t4 %>%
  filter(celltype %in% c("Progenitors" ,"CD8 Exhausted")) %>%
  ggplot(aes(x = survival, y = donor_percentage, color = survival)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +  
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) + 
  facet_wrap(~ celltype, scales = "free_y") +
  stat_compare_means(aes(group = survival), method = "wilcox.test", label.y = 50, label.x = 1.25, size= 3, label="p.format") +  # Mann-Whitney U
  theme_bw() +
  labs(
       x = "Cohort",
       y = "Donor Chimerism") +
  scale_color_manual(values = survival_colors) +
  scale_y_continuous(limits = c(0, 100)) + 
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
    axis.ticks = element_line(color = "black", size = 0.5) 
  ) 


pdf("6.4_donor_chimerism_comparision_prog_exh_MT_only.pdf", width = 8, height = 6)
donor_chimerism_comparision
dev.off()


