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
setwd("~/TP53_ImmuneEscape/5_Souporcell/")

# Load the metadata
seu <- readRDS("~/250128_seurat_annotated_final.rds")

seu@meta.data <- seu@meta.data %>%
  mutate(timepoint = ifelse(timepoint == "pre-transplant", "0", timepoint))

# List of patients
patient_list <- unique(seu$patient_id) %>% sort()

# Initialize a list to store individual patient data
results <- list()

# Initialize a final data frame to store all patients together
final_dataset <- tibble()  # Empty tibble to store merged data

# Loop through each patient
for (patient_id in patient_list) {
  
  # Notifier: Show which patient is being processed
  message(paste("🔄 Processing patient:", patient_id))
  
  patient_meta <- seu@meta.data %>%
    as_tibble(rownames = "barcode") %>%
    filter(patient_id == !!patient_id)
  
  # Load the corresponding Souporcell output file
  souporcell_file <- paste0("clusters/", patient_id, "_clusters.tsv")
  
  # Check if the file exists before proceeding
  if (!file.exists(souporcell_file)) {
    message(paste("⚠️ Warning: Souporcell file not found for", patient_id))
    next  # Skip this patient
  }
  
  df_souporcell <- read_tsv(souporcell_file)
  
  # Wrangle the data frame 
  patient_meta$barcode <- gsub(".*_", "", patient_meta$barcode)  # Extract barcode after last "_"
  
  # Merge metadata with Souporcell data
  df_merged <- patient_meta %>% left_join(df_souporcell, by = "barcode") %>% 
    filter(assignment %in% c("0", "1"))
  
  # ✅ Check if timepoint == 0 sample exists first
  pretransplant_data <- df_merged %>%
    filter(timepoint == 0) %>%
    count(assignment, name = "pretransplant_cells") %>%
    mutate(percent = pretransplant_cells / sum(pretransplant_cells) * 100)
  
  donor_assignment <- NULL
  assignment_source <- "None"
  donor_percentage <- NA_real_
  
  if (nrow(pretransplant_data) > 0) {
    dominant_assignment <- pretransplant_data %>%
      filter(percent > 80) %>%
      pull(assignment)
    
    if (length(dominant_assignment) == 1) {
      donor_assignment <- setdiff(c("0", "1"), dominant_assignment)
      assignment_source <- "Timepoint = 0 sample"
      donor_percentage <- max(pretransplant_data$percent)
    } else {
      message(paste("⚠️ Could not determine dominant genotype from timepoint = 0 sample for", patient_id))
    }
  } else {
    message(paste("⚠️ No pre-transplant sample found for", patient_id))
  }
  
  # Compute total Mid Ery & Late Ery cell count per assignment
  cell_counts <- df_merged %>%
    filter(celltype %in% c("Mid Erythroids", "Late Erythroids")) %>%
    count(assignment, name = "total_ery_cells")
  
  total_ery_cells <- sum(cell_counts$total_ery_cells, na.rm = TRUE)
  
  if (is.null(donor_assignment) && total_ery_cells > 200) {
    if (nrow(cell_counts) > 0) {
      cell_counts <- cell_counts %>%
        mutate(percent = total_ery_cells / sum(total_ery_cells) * 100)
      
      if (max(cell_counts$percent, na.rm = TRUE) > 80) {
        donor_assignment <- cell_counts %>%
          arrange(desc(percent)) %>%
          slice(1) %>%
          pull(assignment)
        assignment_source <- "Mid/Late Erythroid Ratio"
        donor_percentage <- max(cell_counts$percent)
      } else {
        message(paste("⚠️ Mid/Late Erythroid ratio below 80% or insufficient dominance for", patient_id, "- Moving to all cell type ratio"))
      }
    }
  } else {
    message(paste("⚠️ Total Mid/Late Erythroid cells ≤ 200 for", patient_id, "- Skipping this method"))
  }
  
  if (is.null(donor_assignment)) {
    all_cells_count <- df_merged %>%
      count(assignment, name = "total_all_cells")
    
    if (nrow(all_cells_count) > 0) {
      all_cells_count <- all_cells_count %>%
        mutate(percent = total_all_cells / sum(total_all_cells) * 100)
      
      all_cells_donor_assignment <- all_cells_count %>%
        filter(percent > 80) %>%
        pull(assignment)
      
      
      if (length(all_cells_donor_assignment) == 1) {
        donor_assignment <- all_cells_donor_assignment
        assignment_source <- "Overall Cell Type Ratio"
        donor_percentage <- max(all_cells_count$percent)
      } else {
        # 🔥 **Fix: Assign "Unknown" when no dominant donor is found**
        donor_assignment <- "unknown"
        assignment_source <- "unknown"
        donor_percentage <- NA_real_
        message(paste("⚠️ No dominant donor found for", patient_id, "- Assigning as unknown"))
      }
    }
  }
  
  if (!is.null(donor_assignment)) {
    # 🔥 **Fix: Ensure 'Unknown' is correctly assigned in the dataset**
    df_merged <- df_merged %>%
      mutate(assignment = as.character(assignment),
             origin = case_when(
               donor_assignment == "unknown" ~ "unknown",  # Assign unknown cases
               assignment == donor_assignment ~ "donor",
               TRUE ~ "recipient"
             ),
             patient_id = patient_id)
    
    # Store patient data in results list
    results[[patient_id]] <- df_merged
    
    # Append to final dataset
    final_dataset <- bind_rows(final_dataset, df_merged)
    
    message(paste("✅ Assigned donor using", assignment_source, "for", patient_id, "- Donor Genotype:", donor_assignment, "(", round(donor_percentage, 2), "% )"))
  } else {
    message(paste("No clear donor assignment for", patient_id, "- Assigned as Unknown"))
  }
}

# Save final dataset as CSV
if (nrow(final_dataset) > 0) {
  write_csv(final_dataset, "final_dataset.csv")
  message("📁 Final dataset saved as final_dataset.csv")
} else {
  message("⚠️ No valid patient data processed.")
}

# Visualize the results 

# Filter for Erythroid cells 
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
  facet_wrap(~celltype, scales = "free_x") + # Separate panels for Mid & Late Erythroids
  theme_minimal() +
  labs(title = "Genotype Proportions in Remission Samples",
       x = "Sample ID",
       y = "Proportion",
       fill = "Genotype") +
  #scale_fill_manual(values = c("0" = "blue", "1" = "red")) + # Custom colors
  scale_fill_manual(values = c("donor" = "#377eb8", 
                               "recipient" = "#c44e52", 
                               "unknown" = "#b0b0b0")) +
  # Improve X-axis readability & bring back ticks
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), # Rotate 45 degrees
        axis.ticks.x = element_line(), # Bring back x-axis ticks
        strip.text = element_text(size = 10, face = "bold"), # Make facet titles bigger
        panel.spacing = unit(1.5, "lines"), # Increase spacing between facets
        legend.position = "right") +  # Move legend to the right for more space
  
  # Show every sample ID but only the first 3 characters
  scale_x_discrete(labels = function(x) substr(x, 1, 3),
                   expand = c(0.05, 0.05))
p1
# Save as a pdf file 
pdf("5.4_souporcell_results_.pdf", width = 12, height = 8)
p1
dev.off()

#################################

# Subset the corresponding patient from the metadata        
P02 <- subset(x = seu, subset = patient_id=="P02")
        
# Load the souporcell output tsv file that contains the assignments from the souporcell
Df1 <- read_tsv("outputs/clusters/P02_clusters.tsv")

# Wrangle the data frame 
P02 = as.data.frame(P02@meta.data) %>% rownames_to_column("barcode")
P02$barcode <- gsub(".*_", "", P02$barcode)

# Merge 2 data frames  by using the joint column.
df02 <- P02 %>% left_join(Df1, by="barcode")

# Check each sample with their souporcell assignment and decide on the host and donor cells using the pretransplant sample as a guide. 
subset(df02, subset = assignment %in% c("0","1")) %>% 
  tabyl(celltype, assignment,sample_id, show_missing_levels = FALSE) %>%
  adorn_totals("row")%>% 
 adorn_percentages("row") %>%  # Add percentages within each row
  adorn_rounding(digits = 2) %>%  # Round to 2 decimal places
  adorn_pct_formatting()

# Only keep the correct assignments
df02 <- subset(df02, subset = assignment %in% c("0","1")) 

# Rename the assignments, checking with pre-transplant samples which would be host cells
df02$assignment <- gsub("1","host", df02$assignment)
df02$assignment <- gsub("0", "donor", df02$assignment) 

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
