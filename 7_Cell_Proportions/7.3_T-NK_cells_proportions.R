# Nurefsan Sariipek, 231103 updated at 250225
# Check T/NK cell proportions
# Load the libraries
library(ggalluvial)
library(tidyverse)
library(stringr)
library(janitor)
library(ggforce)
library(RColorBrewer)
library(randomcoloR)
library(ggrepel)
library(ggpubr)
library(ggnewscale)

# Empty environment
rm(list=ls())

# For Nurefsan:
setwd("~/TP53_ImmuneEscape/7_Cell_Proportions/")

# Load the saved dataframe that contains the information about souporcell information
# For Nurefsan
final_df <- read_csv("~/final_dataset.csv")

# final_df <- final_df %>%
#   mutate(TP53_status = ifelse(as.numeric(str_extract(patient_id, "\\d+")) %in% c(1:12, 14, 17), "MT", "WT"))
# #save the version with mut info 
# 
# write_csv(final_df, "~/final_dataset.csv")

# Calculate the proportions for T and NK cells
proportions_df <- final_df %>%
  filter(celltype %in% c("CD4 Naïve","CD4 Memory","CD4 Effector Memory","Treg","CD8 Naïve","CD8 Memory",
                         "CD8 Effector","CD8 Exhausted", "γδ T","NK T","CD56 Dim NK", 
                         "CD56 Bright NK","Adaptive NK") &
        timepoint %in% c("3","5","6") & 
        sample_status =="remission" 
          # & TP53_status=="MT"
        # & origin == "donor"
       ) %>% 
  group_by(sample_id, TP53_status) %>% reframe(tabyl(celltype)) %>%
  mutate(percent = percent*100) %>% 
  mutate(celltype = factor(celltype,
                           levels = c("CD4 Naïve","CD4 Memory", "CD4 Effector Memory" ,"Treg","CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Exhausted", "γδ T","NK T","CD56 Dim NK", "CD56 Bright NK", "Adaptive NK")))

# # Calculate the proportions for myeloid compartments
# proportions_df <- final_df %>%
#   filter(celltype %in% c("Progenitors","Early Erythroids","Mid Erythroids" ,"Late Erythroids","pDC","cDC","Pro Monocytes", "Monocytes","Non Classical Monocytes") &
#            timepoint %in% c("3","5","6") & sample_status =="remission" 
#           # origin == "donor"
#   ) %>% 
#   group_by(sample_id, survival) %>% reframe(tabyl(celltype)) %>%
#   mutate(percent = percent*100) %>% 
#   mutate(celltype = factor(celltype,
#                            levels = c("Progenitors","Early Erythroids","Mid Erythroids" ,"Late Erythroids","pDC","cDC","Pro Monocytes", "Monocytes","Non Classical Monocytes")))
              
# Visualize 
p1 <- proportions_df %>%
  mutate(TP53_status = factor(TP53_status, levels = c("MT", "WT"))) %>%
  ggplot(aes(x = TP53_status, y = percent, color = TP53_status)) +
  geom_point(shape = 1, size = 3) +
  coord_cartesian(ylim = c(0,65)) +
  facet_wrap(~ celltype, ncol = 3) +
  theme_bw() +
 # scale_color_manual(values = c("Relapsed"= "tomato1", "Non-relapsed"="royalblue1"))+
  theme(aspect.ratio = 0.75,
        axis.text.x = element_text(face="plain", size=12, color="black"), 
        #axis.title.x = element_text(face="plain", size=20, color="black"),
        axis.text.y = element_text(face="plain", size=12, color="black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size=20, face="plain"),
        strip.text = element_text(size=12, face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14)) +
        theme(panel.grid.minor = element_blank())

# Check the plot
p1
# Save as a pdf
pdf("7.3_T-NK_proportions_WTvsMT.pdf", width = 10, height = 8)
p1
dev.off()


# Statistical tests
print(proportions_df, n = 100)
for (x in unique(proportions_df$celltype)) {
  print(x)
  print(t.test(filter(proportions_df, celltype == x, survival == "Non-relapsed")$percent,
               filter(proportions_df, celltype == x, survival == "Relapsed")$percent,
               var.equal = T)$p.value)
}


# Heatmap
p2 <- proportions_df %>%
  ggplot(aes(x = sample_id, y = celltype, fill = percent)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma") +  # First color scale for percent fill
  ggnewscale::new_scale_fill() +  # Allows another color scale
  geom_tile(aes(x = sample_id, y = -1, fill = TP53_status), height = 0.1) +  # Small tile below for survival
  #scale_fill_manual(values = c("Relapsed"= "tomato1", "Non-relapsed"="royalblue1")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

p2
# Save as a pdf
pdf("7.3_NK-Tcells_WTvsMT_heatmap.pdf", width = 16, height = 8)
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


