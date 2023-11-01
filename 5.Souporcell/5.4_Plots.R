#load the libraries

library(ggalluvial)

#Load the saved dataframe from the previous part

combined_df <- read_csv(file = "~/cohort1-2_souporcell.csv")



#subset and summarise as you like to visualize 

#select only T cells
combined_df_T <- subset(combined_df, subset = celltype %in% c("CD4 Memory","CD8 Effector","CD8 Memory","CD4 Naïve","Treg","CD8 Naïve","CD8 Effector","CD4 Naïve","CD56 Dim NK cells","CD8 Terminally Exhausted","CD4 Memory","NK T cells","γδ T lymphocytes","CD56 Bright NK cells", "CD4 Naïve"))

#select only remission cells
combined_df_T_rem <- subset(combined_df_T, subset = status == "remission")

#select only 6mo remission cells and MNC libraries only to prevent any skewing
combined_df_T_rem6mo <- subset(combined_df_T_rem, subset = id %in% c("P01.1Rem","P01.2Rem","P02.1Rem","P05.1Rem
  "."P06.1Rem","P07.1Rem","P08.1Rem"))

#Summarize for plotting purposes
 meta_summary <- 
   combined_df_T_rem6mo  %>% group_by(celltype, cohort,id,assignment) %>%
   summarize(n_cells = n())


ordered_status = c("pre_transplant","remission","relapse")
meta_summary$status= factor(meta_summary$status, levels = ordered_status)

#Plot1
#change fill to change the flow
ggplot(data = meta_summary,
       aes(axis1 = cohort , axis2 = assignment, axis3 = celltype,
           y = n_cells)) +
  scale_x_discrete(limits = c("cohort","origin","celltype"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = celltype)) +
  geom_stratum(aes(fill=celltype), color="black") +
  geom_text(stat = "stratum", aes(label =paste(after_stat(stratum)))) +
  scale_fill_manual(values = palette) + 
  ggtitle("Post-transplant 3-6 months remission samples")












