# Nurefsan Sariipek, 231103
# Load the libraries
library(ggalluvial)
library(tidyverse)
library(janitor)
library(ggforce)
library(RColorBrewer)
library(randomcoloR)
library(ggrepel)
library(fossil)
library(ggpubr)
library(readxl)
# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"

# Load the saved dataframe that contains the information about souporcell information
# For Nurefsan
combined_df <- read_csv(paste0(my_wd,file = "AnalysisNurefsan/Souporcell/outputs/cohort1-2_souporcell.csv"))

# Subset and summarise as you like to visualize 

# Calculate the proportions
proportions_df <- combined_df %>%
  filter(library_type == "MNC",
         celltype %in% c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector", "CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells"
                         ,"CD56 Dim NK cells", "CD56 Bright NK cells"
         ),
         id %in% c("P01.1Rem","P01.2Rem","P02.1Rem","P04.1Rem", "P05.1Rem","P06.1Rem","P07.1Rem","P08.1Rem"),
         assignment == "donor") %>% 
  group_by(id, cohort) %>% reframe(tabyl(celltype)) %>%
  mutate(percent = percent*100) %>% 
  mutate(celltype = factor(celltype,
                           levels = c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Terminally Exhausted", "γδ T lymphocytes","NK T cells","CD56 Dim NK cells", "CD56 Bright NK cells")),
         cohort = gsub("Remission cohort", "Non-relapsed", gsub("Relapse cohort", "Relapsed", cohort)),
         cohort = factor(cohort, levels = c("Non-relapsed", "Relapsed")))

# Visualize 
p1 <- proportions_df %>%
  ggplot(aes(x = cohort, y = percent, color = cohort)) +
  geom_point(shape = 1, size = 3) +
  coord_cartesian(ylim = c(0,50)) +
  facet_wrap(~ celltype, ncol = 3) +
  theme_bw() +
  scale_color_manual(values = c("Relapsed"= "tomato1", "Non-relapsed"="royalblue1"))+
  stat_compare_means(method = "t.test")+
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
pdf("donor_proportions.pdf", width = 10, height = 8)
p1
dev.off()

############################# Peter's additions ###############################

# Statistical tests
print(proportions_df, n = 100)
for (x in unique(proportions_df$celltype)) {
  print(x)
  print(t.test(filter(proportions_df, celltype == x, cohort == "Non-relapsed")$percent,
               filter(proportions_df, celltype == x, cohort == "Relapsed")$percent,
               var.equal = T)$p.value)
}
# CD56 Bright NK cells are significant (P < 0.05)

# View table with percentages
df1<- proportions_df %>%
  pivot_wider(id_cols = "celltype", names_from = "id", values_from = "percent")

View(df1)

#add a new row for all NK cells
df1 %>%
  filter(celltype %in% c("CD56 Bright NK cells", "CD56 Dim NK cells")) %>%
  summarize(celltype = "NK cells", across(c( P01.1Rem, P01.2Rem, P02.1Rem, P04.1Rem, P05.1Rem, P06.1Rem, P07.1Rem, P08.1Rem), sum)) %>%
  bind_rows(df1, .)

# Heatmap
proportions_df %>%
  ggplot(aes(x = id, y = celltype, fill = percent)) +
  geom_tile()
