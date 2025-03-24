# Nurefsan Sariipek, 250220
# Visualize the souporcell results 

library(tidyverse)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(janitor)

# Empty environment
rm(list=ls())

# For Nurefsan:
setwd("~/TP53_ImmuneEscape/5_Souporcell/")

#Load the saved souporcell result table
final_dataset <- read_csv("~/final_dataset.csv")

# Reorder for visulization
final_dataset$sample_status <- factor(final_dataset$sample_status, levels = c("pre_transplant","remission","relapse"))

# Organize the dataset 
t1 <- final_dataset %>%
  filter(sample_status=="remission") %>%
  filter(celltype %in% c("Mid Erythroids", "Late Erythroids")) %>%
  count(sample_id, celltype, sample_status, origin, name = "count") %>%
  group_by(sample_id, celltype) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Create stacked bar plot
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
  scale_fill_manual(values = c("donor" = "#377eb8", 
                               "recipient" = "#c44e52", 
                               "unknown" = "#b0b0b0")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), 
        axis.ticks.x = element_line(), 
        strip.text = element_text(size = 10, face = "bold"), # Make facet titles bigger
        panel.spacing = unit(1.5, "lines"), # Increase spacing between facets
        legend.position = "right") +  # Move legend to the right for more space
  scale_x_discrete(labels = function(x) substr(x, 1, 3),
                   expand = c(0.05, 0.05)) # Show every sample ID but only the first 3 characters

p1
# Save as a pdf file 
pdf("5.5_souporcell_results_erythroids.pdf", width = 12, height = 8)
p1
dev.off()

# Organize the dataset 
t2 <- final_dataset %>%
  filter(sample_status=="remission") %>%
  count(sample_id, celltype, sample_status, origin, name = "count") %>%
  group_by(sample_id, celltype) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Create stacked bar plot
p2 <-  t2 %>%
  ggplot(aes(x = sample_id, y = proportion, fill = origin)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~celltype, scales = "free_x") + # Separate panels for each celltype
  theme_minimal() +
  labs(title = "Genotype Proportions in Remission Samples",
       x = "Sample ID",
       y = "Proportion",
       fill = "Genotype") +
  scale_fill_manual(values = c("donor" = "#377eb8", 
                               "recipient" = "#c44e52", 
                               "unknown" = "#b0b0b0")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), 
        axis.ticks.x = element_line(), 
        strip.text = element_text(size = 10, face = "bold"), # Make facet titles bigger
        panel.spacing = unit(1.5, "lines"), # Increase spacing between facets
        legend.position = "right") +  # Move legend to the right for more space
  scale_x_discrete(labels = function(x) substr(x, 1, 3),
                   expand = c(0.05, 0.05)) # Show every sample ID but only the first 3 characters
p2
# Save as a pdf file 
pdf("5.5_souporcell_results_all_cells_.pdf", width = 12, height = 8)
p2
dev.off()

t2.5 <- final_dataset %>%
  filter(origin %in% c("donor", "recipient")& sample_status == "remission"& !celltype %in% c("UD1", "UD2", "UD3")) %>%
  count(sample_id, celltype, origin, name = "count") %>%  # Count cells
  group_by(sample_id, celltype) %>%  # Group by sample and cell type
  mutate(total_cells = sum(count), 
         donor_percentage = sum(ifelse(origin == "donor", count, 0)) / total_cells * 100) %>%  
  ungroup() %>%
  select(sample_id, celltype, donor_percentage) %>% distinct()


# Heatmap of proportion donor-derived cells
p2.5 <- ggplot(t2.5, aes(x = celltype, y = sample_id, fill = donor_percentage)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", "pink", "red", "darkred"), 
                       limits = c(0, 100),
                       name="Percentage") + 
  labs(x = "Cell Type", y = "Sample ID",  
       title = "Donor Chimerism by Souporcell in Remission Samples"  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Bigger & bold x-axis labels
    axis.text.y = element_text(size = 14),  # Bigger & bold y-axis labels
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.title.y = element_text(size = 16, face = "bold") ,
    legend.title = element_text(size=14))
p2.5
# Save as a pdf file 
pdf("5.5_souporcell_heatmap_.pdf", width = 14, height = 14)
p2.5
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
  scale_fill_manual(values = c("donor" = "#377eb8", 
                               "recipient" = "#c44e52")) +
  labs(x = "Patient ID", y = "Proportion", fill = "Cell Origin",
       title = "Proportion of Cell Origins (Donor vs. Recipient) in Each Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3

# Save as a pdf file 
pdf("5.5_souporcell_per_patient_.pdf", width = 12, height = 8)
p3
dev.off()


