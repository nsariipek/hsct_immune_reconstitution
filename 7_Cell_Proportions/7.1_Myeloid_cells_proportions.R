# Calculate the myeloid cells proportions
# Load the libraries
library(tidyverse)
library(janitor)
library(ggrepel)
library(ggpubr)

# Empty environment
rm(list=ls())

# For Nurefsan:
setwd("~/TP53_ImmuneEscape/7_Cell_Proportions/")

# Load the saved dataframe that contains the information about souporcell information
# For Nurefsan
final_df <- read_csv("~/final_dataset.csv")

# Calculate the proportions for myeloid compartments
proportions_df <- final_df %>%
  filter(celltype %in% c("Progenitors","Early Erythroids","Mid Erythroids" ,"Late Erythroids","pDC","cDC","Pro Monocytes", "Monocytes","Non Classical Monocytes") &
           timepoint %in% c("3","5","6") & sample_status =="remission"
          # origin == "donor"
  ) %>%
  group_by(sample_id, survival) %>% reframe(tabyl(celltype)) %>%
  mutate(percent = percent*100) %>%
  mutate(celltype = factor(celltype,
                           levels = c("Progenitors","Early Erythroids","Mid Erythroids" ,"Late Erythroids","pDC","cDC","Pro Monocytes", "Monocytes","Non Classical Monocytes")))


# Visualize 
p1 <- proportions_df %>%
  mutate(survival = factor(survival, levels = c("Relapsed", "Non-relapsed"))) %>%
  ggplot(aes(x = survival, y = percent, fill = survival)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8, linewidth = 0.3, color = "black") +
  geom_jitter(shape = 21, size = 1.8, width = 0.2, stroke = 0.2, color = "black", aes(fill = survival)) +
  coord_cartesian(ylim = c(0, 65)) +
  facet_wrap(~ celltype, ncol = 3) +
  scale_fill_manual(values = c("Relapsed" = "#E64B35FF", "Non-relapsed" = "#4DBBD5FF")) +
  labs(y = "% of cells", x = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title.y = element_text(size = 8, face = "plain"),
    strip.text = element_text(size = 8, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "none",
    aspect.ratio = 0.9,
    plot.margin = margin(5, 5, 5, 5)
  )

p1

# Save as a pdf
pdf("7.1_Myeloid_proportions_3-6mo.pdf", width = 10, height = 8)
p1
dev.off()
