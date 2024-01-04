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
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/AnalysisNurefsan/"
# For Peter:

# Load the metadata
seu_diet_merged <- readRDS(paste0(my_wd, "AnalysisNurefsan/RDS files/seu_diet_merged.rds"))

# Subset the patient9 from the metadata        
P12 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("9355_MNC","1013_MNC"))
        
# Load the souporcell output tsv file that contains the assignments from the souporcell
Df1 <- read_tsv(paste0(my_wd, "Souporcell/output/clusters/souporcell_Pt11_clusters.tsv"))


# Wrangle the data frame 
P12 = as.data.frame(P12@meta.data) %>% rownames_to_column("barcode")
P12$barcode <- gsub("_.*", "", P12$barcode)

# Merge 2 data frames  by using the joint column.
df11 <- P11 %>% left_join(Df1, by="barcode")

# Check each sample with their souporcell assignment and decide on the host and donor cells using the pretransplant sample as a guide. 
df11 %>%
  tabyl(celltype, assignment,id , show_missing_levels = FALSE) %>%
  adorn_totals("row")

# Only keep the correct assignments
df11 <- subset(df11, subset = assignment %in% c("0","1")) 

# Rename the assignments, checking with pre-transplant samples which would be host cells
df11$assignment <- gsub("0","host", df11$assignment)
df11$assignment <- gsub("1", "donor", df11$assignment) 


#### Repeat this for each patient repeatedly, ask Peter or Ksenia how to make this in a loop

#combine each cohorts seperately 


#Cohort 1= P01 and P02
combined_df1 <- bind_rows(df1,df2)
write_csv(combined_df1, "~/cohort1_souporcell.csv")
combined_df1 <- read_csv(file = "~/cohort1_souporcell.csv")


#Cohort 2 = P05,P06,P07

combined_df2 <- bind_rows(df5,df6,df8)
write_csv(combined_df2, "~/cohort2_souporcell.csv")
combined_df2 <- read_csv(file = "~/cohort2_souporcell.csv")

write_csv(df3,paste0(my_wd, "/Souporcell/output/cohort3_souporcell.csv"))


#Cohort 3= P09, P10, P11, P12

combined_df3 <- bind_rows(df9,df10,df11,df12)




#For patient03 and patient07 we did not have pre-transplant samples that's why we could not run them on the souporcell and since they had >95% chimerism it was safe to assume they were all donor origin at this time point.

#cohort 1- remission cohort
  
#load the previously saved combined data frame which includes patient01 and patient 02 and select the only revelent columns

combined_df1 <- read_csv(file = "~/cohort1_souporcell.csv")

df1 <- 
combined_df1 %>%
dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status.x, assignment, cohort)

#rename the status.x object since souporcell changed to .x and .y
df1$status <- df1$status.x
df1 <- df1 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment, cohort)
  
View(P3)

#wrangle the df
P3 = as.data.frame(P3@meta.data) %>% rownames_to_column("barcode")
#add a column that will match with souporcell output
P3$assignment <- c("donor")

#select the same columns as the other data frame that you'll end up combine 
df3 <- P3 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment,cohort)

#add P3 to cohort 1 dataframe
combined_df1 <- bind_rows(df1,df3)

#now combined_df1 has P01, P02 and P03 information


#Cohort 2- relapse cohort 

#load the previously saved combined data frame which includes patient01 and patient 02 and select the only revelent columns

combined_df2 <- read_csv(file = "~/cohort2_souporcell.csv")

df3 <- 
combined_df3 %>%
dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status.x, assignment, cohort)

#rename the status.x object since souporcell changed to .x and .y
df3$status <- df3$status.x
df3 <- df3 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment, cohort)
  
View(P7)

#wrangle the df
P7 = as.data.frame(P7@meta.data) %>% rownames_to_column("barcode")
#add a column that will match with souporcell output
P7$assignment <- c("donor")

#select the same columns as the other data frame that you'll end up combine 
df7 <- P7 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment,cohort)

#add P7 to cohort 2 dataframe
combined_df2 <- bind_rows(df2,df7)
#now combined_df2 has P05, P06, P07 and P08 information


#combine all of them in a one dataframe and save the new dataframe to use in the future

combined_df<- bind_rows(combined_df1,combined_df2)


#rename the cohorts name so it'll be easier to understand by the audience
combined_df$cohort <- gsub("cohort1", "Remission cohort", combined_df$cohort) 
combined_df$cohort <- gsub("cohort2", "Relapse cohort", combined_df$cohort) 


write_csv(combined_df, "~/cohort1-2_souporcell.csv")
combined_df <- read_csv(file = "~/cohort1-2_souporcell.csv")

#select the remission samples between 3-6 months post-HSCT and T cells only since we are interested in them mostly.

#select only T cells
combined_df_T <- subset(combined_df, subset = celltype %in% c("CD4 Memory","CD8 Effector","CD8 Memory","CD4 Naïve","Treg","CD8 Naïve","CD8 Effector","CD4 Naïve","CD56 Dim NK cells","CD8 Terminally Exhausted","CD4 Memory","NK T cells","γδ T lymphocytes","CD56 Bright NK cells", "CD4 Naïve"))

#select only remission cells
combined_df_T_rem <- subset(combined_df_T, subset = status == "remission")

#select only 6mo remission cells and MNC libraries only to prevent any skewing
combined_df_T_rem6mo <- subset(combined_df_T_rem, subset = id %in% c("P01.1Rem","P01.2Rem","P02.1Rem","P05.1Rem","P06.1Rem","P07.1Rem","P08.1Rem"))

#select the donor cells only for the statistical analysis purposes

combined_df_T_rem6mo_donor <- subset(combined_df_T, subset = assignment == "donor")

#Export as a CSV to do statistical anaylsis

t1 <- combined_df_T_rem6mo_donor %>% 
  group_by(celltype, cohort,id) %>%
  summarize(n = n()) %>%
  group_by(id,cohort) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n)

head(t1)
 
write_csv(t1, "~/combined_df_T_rem_donor.csv")


            
 
