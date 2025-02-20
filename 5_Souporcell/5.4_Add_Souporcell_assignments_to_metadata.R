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



# Save final dataset as CSV, saved this to auxillary file on the dropbox since it is too big to snyc at github
if (nrow(final_dataset) > 0) {
  write_csv(final_dataset, "final_dataset.csv")
  message("📁 Final dataset saved as final_dataset.csv")
} else {
  message("⚠️ No valid patient data processed.")
}

# Visualize the results 
final_dataset$sample_status <- factor(final_dataset$sample_status, levels = c("pre_transplant","remission","relapse"))

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


t2 <- final_dataset %>%
  filter(origin %in% c("donor", "recipient")) %>%
  count(patient_id,origin,sample_status, name = "count") %>%
  group_by(patient_id, sample_status) %>%  
  mutate(proportion = count / sum(count)) %>%
  ungroup()


p3 <- ggplot(t2, aes(x = sample_status, y = proportion, fill = origin)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ patient_id) +
  scale_fill_manual(values = c("donor" = "#377eb8", 
                               "recipient" = "#c44e52")) +
  labs(x = "Patient ID", y = "Proportion", fill = "Cell Origin",
       title = "Proportion of Cell Origins (Donor vs. Recipient) in Each Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save as a pdf file 
pdf("5.4_souporcell_per_patient_.pdf", width = 12, height = 8)
p3
dev.off()



