# Nurefsan Sariipek, 231001
# Distuinguishing donor and host cells using Souporcell

library(tidyverse)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(janitor)


# Load the metadata
seu_diet_merged <- readRDS("~/seu_diet_merged.rds")

# Subset the patient9 from the metadata        
P9 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("1677_MNC","1677_CD3","1732_MNC","1811_MNC","1811_CD3"))
        
# Load the souporcell output tsv file that contains the assignments from the souporcell
Df1 <- read_tsv("~/souporcell_Pt9_clusters.tsv")

# Wrangle the data frame 
P9 = as.data.frame(P9@meta.data) %>% rownames_to_column("barcode")
P9$barcode <- gsub("_.*", "", P9$barcode)

# Merge 2 data frames  by using the joint column.
df9 <- P9 %>% left_join(Df1, by="barcode")

# Check each sample with their souporcell assignment and decide on the host and donor cells using the pretransplant sample as a guide. 
df9 %>%
  tabyl(celltype, assignment,id , show_missing_levels = FALSE) %>%
  adorn_totals("row")

# Only keep the correct assignments
df9 <- subset(df9, subset = assignment %in% c("0","1")) 

# Rename the assignments
df9$assignment <- gsub("1", "donor", df9$assignment)
df9$assignment <- gsub("0", "host", df9$assignment) 


combined_df1 <- bind_rows(df1,df2)
combined_df1_rem <- subset(combined_df1, subset = status.x %in% c("remission"))

combined_df2 <- bind_rows(df5,df6,df8)
combined_df2_rem <- subset(combined_df2, subset = status.x %in% c("remission"))


combined_df1_rem_6mo <- subset(combined_df1_rem, subset = id %in% c("P01.1Rem","P01.1RemT","P02.1Rem","P01.2Rem"))
combined_df2_rem_6mo <- subset(combined_df2_rem, subset = id %in% c("P05.1Rem","P06.1Rem","P08.1Rem","P08.1RemT"))


combined_df1_rem_6mo_donor<- subset(combined_df1_rem_6mo, subset = assignment == "donor")
combined_df2_rem_6mo_donor<- subset(combined_df2_rem_6mo, subset = assignment == "donor")


combined_df1_rem_6mo_host<- subset(combined_df1_rem_6mo, subset = assignment == "host")
combined_df2_rem_6mo_host<- subset(combined_df2_rem_6mo, subset = assignment == "host")


combined_donor6mo <- bind_rows(combined_df1_rem_6mo_donor, combined_df2_rem_6mo_donor)

combined_host6mo <- bind_rows(combined_df1_rem_6mo_host, combined_df2_rem_6mo_host)

combined_df1_donor <-  subset(combined_df1_rem, subset = assignment == "donor")
combined_df2_donor <-  subset(combined_df2_rem, subset = assignment == "donor")

combined_donor <- bind_rows(combined_df1_rem_6mo_donor, combined_df2_rem_6mo_donor)


#to export the table 
t3 <-
  combined_donor %>%
  group_by( celltype, id) %>% 
  summarize(n=n()) 

write_csv(t3, "~/donor.csv")

              
# Visualize the results in different ways
 

# Visualize the souporcell outputs at which time point 
df4 %>% 
  select(assignment, id) %>% 
  group_by(assignment, id) %>% 
  summarize(n=n()) %>% 
  ggplot(aes(x= assignment, y=n, fill=id, label= n)) + 
  geom_bar(stat = "identity") + 
  geom_text(size = 3, position = position_stack(vjust = 0.5),colour = "white")+
  theme_pubr(base_size = 15) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5), legend.position = "right")


combined_donor  %>% 
     select(celltype,cohort) %>% 
     group_by(celltype,cohort) %>% 
     summarize(n=n()) %>% 
     ggplot(aes(x= cohort, y=n, fill=celltype)) + 
     geom_bar(stat = "identity", position = "fill") + 
     theme_pubr(base_size = 15) +
     scale_fill_manual(values = palette) +
     theme(axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10),legend.position = "right") 



combined_donor %>% 
   select(id, celltype) %>% 
   group_by(id, celltype) %>% 
   summarize(n=n()) %>% 
   ggplot(aes(x= celltype, y=n, fill=id, label= n)) + 
   geom_bar(stat = "identity") + 
   geom_text(size = 3, position = position_stack(vjust = 0.5),colour = "white")+
   theme_pubr(base_size = 15) +
   scale_fill_manual(values = palette)+ 
  theme(axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10), legend.position = "right")
