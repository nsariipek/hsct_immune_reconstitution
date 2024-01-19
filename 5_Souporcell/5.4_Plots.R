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


#Load the saved dataframe that contains the information from the previous part
# For Nurefsan
combined_df <- read_csv(file = "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/AnalysisNurefsan/Souporcell/outputs/cohort1-2_souporcell.csv")
# For Peter
combined_df <- read_csv(file = "~/DropboxMGB/Projects/ImmuneEscapeTP53/AnalysisNurefsan/Souporcell/outputs/cohort1-2_souporcell.csv")

#subset and summarise as you like to visualize 


#select only T and NK cells
combined_df_T <- subset(combined_df, subset = celltype %in% c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Dim NK cells","CD56 Bright NK cells")) 

#select only remission cells
combined_df_T_rem <- subset(combined_df_T, subset = status == "remission")

#select only 6mo remission cells and MNC libraries only to prevent any skewing
combined_df_T_rem6mo <- subset(combined_df_T_rem, subset = id %in% c("P01.1Rem","P01.2Rem","P02.1Rem","P04.1Rem", "P05.1Rem","P06.1Rem","P07.1Rem","P08.1Rem"))


#Plot 1, Post-transplant 3-6 months remission samples

#Summarize the dataframe
meta_summary <- 
  combined_df_T_rem6mo  %>% group_by(celltype,cohort,id,assignment) %>%
  mutate(cohort = gsub("Remission cohort", "Non-relapsed", gsub("Relapse cohort", "Relapsed", cohort)),
         cohort = factor(cohort, levels = c("Non-relapsed", "Relapsed")))  %>%
  summarize(n_cells = n())
#add a color palette 
palette <- distinctColorPalette(k = 34)


#change fill to change the flow
ggplot(data = meta_summary,
       aes(axis1 = cohort , axis2 = assignment, axis3 = celltype,
           y = n_cells)) +
  scale_x_discrete(limits = c("Cohort","Origin","Celltype"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = celltype)) +
  geom_stratum(aes(fill=celltype), color="black") +
  geom_text(stat = "stratum", aes(label =paste(after_stat(stratum))), size =5) +
  scale_fill_manual(values = palette) + 
  theme(axis.text.x = element_text(face="bold", size=16, color="black"), 
        axis.title.x = element_text(face="bold", size=16, color="black"),
        axis.text.y = element_text(face="plain", size=12, color="black"),
        axis.title.y = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        plot.title = element_text(size=20, face="plain"),
        legend.position = "none") +
  ggtitle("Post-transplant 3-6 months remission samples")

  

# Peter's quick additions
# Bar/dot plots
proportions_df <- combined_df %>% filter(library_type == "MNC",
                       celltype %in% c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector",
                                       "CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Dim NK cells",
                                       "CD56 Bright NK cells"),
                       id %in% c("P01.1Rem","P01.2Rem","P02.1Rem","P04.1Rem", "P05.1Rem","P06.1Rem","P07.1Rem","P08.1Rem"),
                       status == "remission", # redundant
                       assignment == "donor") %>% 
  group_by(id, cohort) %>% reframe(tabyl(celltype)) %>%
  mutate(percent = percent*100)

proportions_df <- proportions_df %>% mutate(celltype = factor(celltype,
  levels = c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Terminally Exhausted",
             "γδ T lymphocytes","NK T cells","CD56 Dim NK cells","CD56 Bright NK cells")),
  cohort = gsub("Remission cohort", "Non-relapsed", gsub("Relapse cohort", "Relapsed", cohort)),
  cohort = factor(cohort, levels = c("Non-relapsed", "Relapsed")))

#additions by Nurefsan , adding the p values to the plot and making it bigger 
proportions_df %>%
  ggplot(aes(x = cohort, y = percent, color = cohort)) +
  geom_point(shape = 1, size = 2) +
  coord_cartesian(ylim = c(0,50)) +
  facet_wrap(~ celltype) +
  theme_bw() +
  scale_color_manual(values = c("Relapsed"="red", "Non-relapsed"="green"))+
  stat_compare_means(method = "t.test")+
  theme(axis.text.x = element_text(face="plain", size=16, color="black"), 
        axis.title.x = element_text(face="bold", size=20, color="black"),
        axis.text.y = element_text(face="plain", size=12, color="black"),
        axis.title.y = element_blank(),
        plot.title = element_text(size=20, face="plain"),
        strip.text = element_text(size=18, face="plain"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14)) +
       theme(panel.grid.minor = element_blank())

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


#Plot 3, see how accurate souporcell results comparing with chimerism information

df <- read_excel("/Dropbox/ImmuneEscapeTP53/AnalysisNurefsan/Souporcell/output/similarity.xlsx")
View(df)


ids = c("P01.1Rem",  "P01.1RemT", "P01.2Rem" , "P03.1Rem" , "P04.1Rem",  "P04.1RemT", "P05.1Rem" , "P05.Rel" ,  "P05.RelT" ,"P06.1Rem",  "P07.1Rem",  "P08.1Rem" , "P08.2Rem" ,"P09.1Rem", "P09.Rel", "P09.RelT", "P10.1Rem", "P10.Rel", "P11.1Rem" , "P11.1RemT", "P12.1Rem" )
df$id = factor(df$id, levels= ids)


df %>% 
  ggplot(aes(x = Souporcell, y = Chimerism, colour = id)) +
  geom_point(size=3.5, position = "jitter")+
  #scale_color_distiller(palette = "Spectral")+
  #geom_text(hjust=0, vjust=0)+ also add label= id inside of aes 
  #theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  coord_cartesian(expand = FALSE, xlim = c(0, 1.1), ylim = c(0, 1.1))+
  theme(axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.title = element_text(color = "grey20", size = 16, hjust = .5, vjust = .5, face = "plain"))+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  #ggtitle("Correlation of Souporcell Results with Clinical Chimerism Data")+
  theme(aspect.ratio=1)



