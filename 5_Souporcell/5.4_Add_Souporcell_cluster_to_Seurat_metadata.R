# Nurefsan Sariipek, 231001
# Distuinguishing donor and host cells using Souporcell

library(tidyverse)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(janitor)

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"

# Load the metadata
seu_diet_merged <- readRDS(paste0(my_wd, "/RDS files/seu_diet_merged.rds"))

##Note: Nurefsan did these for each patient since she is really bad with loops, in future turn this into a loop, lines down below are representative for just one patient##

# Subset the corresponding patient from the metadata        
P11 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("5641_MNC","6244_MNC","6244_CD3"))
        
# Load the souporcell output tsv file that contains the assignments from the souporcell
Df1 <- read_tsv(paste0(my_wd, "/TP53_ImmuneEscape/5_Souporcell/outputs/clusters/souporcell_Pt11_clusters.tsv"))

# Wrangle the data frame 
P11 = as.data.frame(P11@meta.data) %>% rownames_to_column("barcode")
P11$barcode <- gsub("_.*", "", P11$barcode)

# Merge 2 data frames  by using the joint column.
df11 <- P11 %>% left_join(Df1, by="barcode")

# Check each sample with their souporcell assignment and decide on the host and donor cells using the pretransplant sample as a guide. 
df10 %>%
  tabyl(celltype, assignment,id , show_missing_levels = FALSE) %>%
  adorn_totals("row")

# Only keep the correct assignments
df10 <- subset(df10, subset = assignment %in% c("0","1")) 

# Rename the assignments, checking with pre-transplant samples which would be host cells
df10$assignment <- gsub("1","host", df10$assignment)
df10$assignment <- gsub("0", "donor", df10$assignment) 

##### Combine each cohorts separately ##########

# Cohort 1= P01 and P02
combined_df1 <- bind_rows(df1,df2)
write_csv(combined_df1, "~/cohort1_souporcell.csv")

# Cohort 2 = P05,P06,P07

combined_df2 <- bind_rows(df5,df6,df8)
write_csv(combined_df2, "~/cohort2_souporcell.csv")

# Cohort 3= P09, P10, P11, P12
#Load 
combined3 <- read_csv(paste0(my_wd, "/Trash/cohort3_souporcell.csv"))

combined3 <- combined3[combined3$patient_identity != 'pt10', ]

df10$status <- df10$status.x
df10 <- df10%>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, barcode, id, status, assignment, cohort) %>%
  rename(cell = barcode)

combined_df3 <- bind_rows(df10,combined3)
# Somehow I must have removed the cell identifier for Patient 10 when I re-ran the souporcell which caused a lot of problems in the following analysis so I manually added the same identifier from seu object, which is insanely problematic I know. I'll re-run pt 10 and this time use barcodes with identifier from original dataset.
# combined_df3$cell <- ifelse(combined_df3$orig.ident == "1195_MNC", paste0(combined_df3$cell, "_2"), combined_df3$cell)
write_csv(combined_df3, "~/cohort3_souporcell.csv")

#For patient03 and patient07 we did not have pre-transplant samples that's why we could not run them on the souporcell and since they had >95% chimerism it was safe to assume they were all donor origin at this time point, so we added that manually over here#

#Cohort 1- remission cohort
tb1 <- 
combined_df1 %>%
dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status.x, assignment, cohort)

# Rename the status.x object since souporcell changed to .x and .y
tb1$status <- tb1$status.x
tb1 <- tb1 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment, cohort)


# Subset the corresponding patient from the metadata        
P3 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("",""))  
View(P3)

# Wrangle the df
P3 = as.data.frame(P3@meta.data) %>% rownames_to_column("barcode")
# Add a column that will match with souporcell output
P3$assignment <- c("donor")

# Select the same columns as the other data frame that you'll end up combine 
df3 <- P3 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment,cohort)

# Add P3 to cohort 1 dataframe
combined_df1 <- bind_rows(tb1,df3)
write_csv(combined_df1, "~/cohort1_souporcell.csv")
# Now combined_df1 has P01, P02 and P03 information

# Cohort 2- relapse cohort 

# Load the previously saved combined data frame which includes patient01 and patient 02 and select the only relevant columns
combined_df2 <- read_csv(file = "~/cohort2_souporcell.csv")

tb2 <- 
combined_df2 %>%
dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status.x, assignment, cohort)

# Rename the status.x object since souporcell changed to .x and .y
tb2$status <- tb2$status.x
tb2 <- tb2 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment, cohort)


# Subset the corresponding patient from the metadata        
P7 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("",""))  
View(P7)

# Wrangle the df
P7 = as.data.frame(P7@meta.data) %>% rownames_to_column("barcode")
# Add a column that will match with souporcell output
P7$assignment <- c("donor")

# Select the same columns as the other data frame that you'll end up combine 
df7 <- P7 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment,cohort)

# Add P7 to cohort 2 dataframe
combined_df2 <- bind_rows(tb22,df7)
# Now combined_df2 has P05, P06, P07 and P08 information

# Combine all of them in a one dataframe and save the new dataframe to use in the future
combined_df<- bind_rows(combined_df1,combined_df2)

# Rename the cohorts name so it'll be easier to understand by the audience
combined_df$cohort <- gsub("cohort1", "Remission cohort", combined_df$cohort) 
combined_df$cohort <- gsub("cohort2", "Relapse cohort", combined_df$cohort) 

write_csv(combined_df, "~/cohort1-2_souporcell.csv")

## From now on only use this combined df for the upcoming analysis ##
combined_df <- read_csv(paste0(my_wd, file = "TP53_ImmuneEscape/5_Souporcell/results/cohort1-2_souporcell.csv"))

# # Select the remission samples between 3-6 months post-HSCT and T cells only since we are interested in them mostly.
# # Select only T cells
# combined_df_T <- subset(combined_df, subset = celltype %in% c("CD4 Memory","CD8 Effector","CD8 Memory","CD4 Naïve","Treg","CD8 Naïve","CD8 Effector","CD4 Naïve","CD56 Dim NK cells","CD8 Terminally Exhausted","CD4 Memory","NK T cells","γδ T lymphocytes","CD56 Bright NK cells", "CD4 Naïve"))
# # Select only remission cells
# combined_df_T_rem <- subset(combined_df_T, subset = status == "remission")
# 
# # Select only 6mo remission cells and MNC libraries only to prevent any skewing
# combined_df_T_rem6mo <- subset(combined_df_T_rem, subset = id %in% c("P01.1Rem","P01.2Rem","P02.1Rem","P05.1Rem","P06.1Rem","P07.1Rem","P08.1Rem"))
# 
# # Select the donor cells only for the statistical analysis purposes
# combined_df_T_rem6mo_donor <- subset(combined_df_T, subset = assignment == "donor")
# 
# #Export as a CSV to do statistical anaylsis
# t1 <- combined_df_T_rem6mo_donor %>% 
#   group_by(celltype, cohort,id) %>%
#   summarize(n = n()) %>%
#   group_by(id,cohort) %>%
#   mutate(freq = n/sum(n)) %>%
#   dplyr::select(-n)
# 
# head(t1)
#  
# write_csv(t1, "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Souporcell/results/combined_df_T_rem_donor.csv")
