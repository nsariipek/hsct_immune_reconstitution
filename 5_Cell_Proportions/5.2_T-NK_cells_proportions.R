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

# Load the saved dataframe that contains the information about souporcell information
# For Nurefsan
final_df <- read_csv("~/final_dataset.csv")
final_df$celltype <- gsub("γδ T", "delta-gamma T", final_df$celltype)


# final_df <- final_df %>%
#   mutate(TP53_status = ifelse(as.numeric(str_extract(patient_id, "\\d+")) %in% c(1:12, 14, 17), "MT", "WT"))
# #save the version with mut info 
# 
# write_csv(final_df, "~/final_dataset.csv")

# Calculate the proportions for T and NK cells
proportions_df <- final_df %>%
  filter(celltype %in% c("CD4 Naïve","CD4 Memory","CD4 Effector Memory","Treg","CD8 Naïve","CD8 Memory",
                         "CD8 Effector","CD8 Exhausted", "delta-gamma T") &
        timepoint %in% c("3","5","6") & 
        sample_status =="remission"  
          &  TP53_status=="MT"
        #  & origin == "donor"
       ) %>% 
  group_by(sample_id, survival) %>% reframe(tabyl(celltype)) %>%
  mutate(percent = percent*100) %>% 
  mutate(celltype = factor(celltype,
                           levels = c("CD4 Naïve","CD4 Memory", "CD4 Effector Memory" ,"Treg","CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Exhausted", "delta-gamma T"))) %>%
  mutate(survival = factor(survival, levels = c("Relapsed", "Non-relapsed")))


# Plot
p1 <- ggplot(proportions_df, aes(x = survival, y = percent, fill = survival)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9, linewidth = 0.5, color = "black") +  # Thicker lines
  geom_jitter(shape = 21, size = 1.8, width = 0.15, stroke = 0.2, color = "black", aes(fill = survival)) +
  coord_cartesian(ylim = c(0, 35)) +
  facet_wrap(~ celltype, ncol = 5) +  # ← 4 per row
  scale_fill_manual(values = c("Relapsed" = "#E64B35FF", "Non-relapsed" = "#4DBBD5FF")) +
  labs(y = "% of cells", x = NULL) +
  stat_compare_means(
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
    aspect.ratio = 2,  # ← Keeps it narrow/tall
    plot.margin = margin(4, 4, 4, 4)
  )

# Check the plot
p1
# Save as a pdf
pdf("5.2_T_proportions_TP53_MT_signif.pdf", width = 8, height = 6)
p1
dev.off()

#############################################

# Heatmap
p2 <- proportions_df %>%
  ggplot(aes(x = sample_id, y = celltype, fill = percent)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma") +  # First color scale for percent fill
  ggnewscale::new_scale_fill() +  # Allows another color scale
  geom_tile(aes(x = sample_id, y = -1, fill = survival), height = 0.1) +  # Small tile below for survival
  #scale_fill_manual(values = c("Relapsed"= "tomato1", "Non-relapsed"="royalblue1")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

p2
# Save as a pdf
pdf("7.3_Tcells_MT_heatmap.pdf", width = 16, height = 8)
p2
dev.off()


###############################
# View table with percentages
df1<- proportions_df %>%
  pivot_wider(id_cols = "celltype", names_from = "sample_id", values_from = "percent")

View(df1)

#add a new row for all NK cells
df1 %>%
  filter(celltype %in% c("CD56 Bright NK", "CD56 Dim NK")) %>%
  summarize(celltype = "NK cells", 
            across(c(P01_Rem1, P01_Rem2, P04_Rem, P05_Rem, P06_Rem, P07_Rem,
                     P08_Rem1, P13_Rem, P14_Rem, P15_Rem, P16_Rem, P17_Rem,
                     P18_Rem, P19_Rem, P20_Rem, P21_Rem, P22_Rem, P23_Rem,
                     P24_Rem, P25_Rem, P26_Rem, P27_Rem, P28_Rem, P29_Rem,
                     P30_Rem, P31_Rem, P32_Rem, P33_Rem), sum, na.rm = TRUE)) %>%
  bind_rows(df1, .)


