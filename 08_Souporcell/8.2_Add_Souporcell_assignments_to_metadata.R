# Nurefsan Sariipek and Peter van Galen, 250518
# Add Souporcell results to Seurat metadata and save as gzipped csv file

library(tidyverse)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(janitor)

# Empty environment
rm(list = ls())

# For Nurefsan:
setwd("~/hsct_immune_reconstitution/08_Souporcell/")
# For Peter:
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/08_Souporcell/")

# Load Seurat data
seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")

# List of patients
patient_list <- unique(seu$patient_id) %>% sort()

# Initialize a list to store individual patient data
results <- list()

# Initialize a final data frame to store all patients together
souporcell_assignments <- tibble() # Empty tibble to store merged data

# Loop through each patient
for (patient_id in patient_list) {
  # patient_id <- "P20"
  # Show which patient is being processed
  message(paste("▶️ Processing:", patient_id))

  patient_meta <- seu@meta.data %>%
    as_tibble(rownames = "barcode") %>%
    filter(patient_id == !!patient_id)

  # Load the corresponding Souporcell output file
  souporcell_file <- paste0("clusters/", patient_id, "_clusters.tsv")

  # Check if the file exists before proceeding
  if (!file.exists(souporcell_file)) {
    message(paste("⚠️ Skipping due to Souporcell file not found:", patient_id))
    next
  }

  df_souporcell <- read_tsv(
    souporcell_file,
    col_types = cols(assignment = col_character()),
    show_col_types = FALSE
  )

  # Extract barcode after last "_"
  patient_meta$barcode <- gsub(".*_", "", patient_meta$barcode)

  # Merge metadata with Souporcell data
  df_merged <- patient_meta %>%
    left_join(df_souporcell, by = "barcode") %>%
    filter(assignment %in% c("0", "1"))

  ########## Rule-based classification #########

  # Check if timepoint == 0 sample exists first
  pretransplant_data <- df_merged %>%
    filter(timepoint == 0) %>%
    count(assignment, name = "pretransplant_cells") %>%
    mutate(percent = pretransplant_cells / sum(pretransplant_cells) * 100)

  donor_assignment <- NULL
  donor_percentage <- NA_real_

  if (nrow(pretransplant_data) > 0) {
    dominant_assignment <- pretransplant_data %>%
      filter(percent > 80) %>%
      pull(assignment)

    if (length(dominant_assignment) == 1) {
      donor_assignment <- setdiff(c("0", "1"), dominant_assignment)
      assignment_source <- "pre-transplant sample"
      donor_percentage <- max(pretransplant_data$percent)
    } else {
      message(paste(
        "⚠️ Could not determine dominant genotype from pre-transplant sample for",
        patient_id
      ))
    }
  } else {
    message(paste("⚠️ No pre-transplant sample found for", patient_id))
  }

  if (is.null(donor_assignment)) {
    all_cells_count <- df_merged %>%
      count(assignment, name = "total_all_cells")

    if (nrow(all_cells_count) > 0) {
      all_cells_count <- all_cells_count %>%
        mutate(percent = total_all_cells / sum(total_all_cells) * 100)

      all_cells_donor_assignment <- all_cells_count %>%
        filter(percent > 75) %>%
        pull(assignment)

      if (length(all_cells_donor_assignment) == 1) {
        donor_assignment <- all_cells_donor_assignment
        assignment_source <- "overall cell type ratio"
        donor_percentage <- max(all_cells_count$percent)
      } else {
        max_percentage <- max(all_cells_count$percent, na.rm = TRUE)
        donor_assignment <- "unknown"
        assignment_source <- "unknown"
        donor_percentage <- NA_real_
        message(paste0(
          "⚠️ No dominant donor found for ",
          patient_id,
          ". Highest genotype percentage: ",
          round(max_percentage, 2),
          "%"
        ))
      }
    }
  }

  df_merged <- df_merged %>%
    mutate(
      assignment = as.character(assignment),
      origin = case_when(
        donor_assignment == "unknown" ~ "unknown",
        assignment == donor_assignment ~ "donor",
        TRUE ~ "recipient"
      ),
      patient_id = patient_id
    )

  # Store patient data in results list
  results[[patient_id]] <- df_merged

  # Append to final dataset
  souporcell_assignments <- bind_rows(souporcell_assignments, df_merged)

  # At the end, replace the current message with:
  if (donor_assignment == "unknown") {
    message(paste0(
      "⚠️ Could not determine donor assignment for ",
      patient_id,
      ". All cells marked as unknown"
    ))
  } else {
    message(paste0(
      "✅ Assigned donor using ",
      assignment_source,
      " for ",
      patient_id,
      ". Donor genotype: ",
      donor_assignment,
      " (",
      round(donor_percentage, 2),
      "%)"
    ))
  }
}

# Save final dataset as a csv file (Gzipped to save space)
souporcell_assignments$cell <- paste0(
  souporcell_assignments$orig.ident,
  "_",
  souporcell_assignments$barcode
)
souporcell_assignments_select <- souporcell_assignments %>%
  select(
    cell,
    cohort,
    sample_status,
    patient_id,
    sample_id,
    TP53_status,
    timepoint,
    celltype,
    assignment,
    origin
  )
write_csv(souporcell_assignments_select, "8.2_Souporcell_assignments.csv.gz")
