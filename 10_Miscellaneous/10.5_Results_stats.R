# Peter van Galen, 250705
# Extract metrics for the paper results section

# Load libraries
library(Seurat)
library(tidyverse)
library(janitor)
library(readxl)

# Set working directory
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/10_Miscellaneous")

# Delete environment variables
rm(list = ls())

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")


# RELAPSE PREDICTION TIME -----------------------------------------------------

# How long after transplant were patients in the relapse cohort diagnosed with relapse?
df <- read_excel("10.2_Timepoints.xlsx")
relapse_after_days_df <- df %>%
  filter(Cohort == "Relapse", Sample_type == "Relapse") %>%
  group_by(Patient_id) %>%
  arrange(Patient_id, Timepoint_days) %>%
  slice_head(n = 1)
summary(as.numeric(relapse_after_days_df$Timepoint_days) / 30.44)

# Subset for six patients with >5% persistent HSPCs ~3 months after transplant (see merged_prog_props_tib from 8.3_Souporcell_plots.R)
subset_relapse_after_days <- relapse_after_days_df %>%
  filter(Patient_id %in% c("P20", "P21", "P22", "P24", "P26", "P27")) %>%
  pull(Timepoint_days) %>%
  as.numeric()
summary(subset_relapse_after_days / 30.44)
# Either way, the median is 5.5 months

# --> "All 10 patients with <5% recipient HSPCs ~3 months after transplant remained in remission for more than three years (remission cohort and patient 23), whereas all six patients with >5% developed relapse within 18 months (median 5.5, range 4.4-16.5 months)."

# Here's another way to get at the same conclusion (more similar to 8.4_Plots.R)
merged_prog_counts_tib <- as_tibble(seu@meta.data) %>%
  filter(
    sample_status == "remission",
    timepoint %in% c(3, 5, 6),
    souporcell_origin %in% c("donor", "recipient"),
    celltype %in%
      c("HSC MPP", "MEP", "LMPP", "Cycling Progenitor", "Early GMP")
  ) %>%
  dplyr::count(patient_id, cohort, souporcell_origin) %>%
  pivot_wider(
    names_from = souporcell_origin,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(donor_percentage = donor / (donor + recipient) * 100) %>%
  arrange(donor_percentage) %>%
  print(n = 30)


# PERSISTENT RECIPIENT PROGENITORS --------------------------------------------

as_tibble(seu@meta.data) %>%
  filter(
    sample_status == "remission",
    timepoint %in% c(3, 5, 6),
    patient_id %in% c("P20", "P21", "P22", "P24", "P26", "P27"),
    !is.na(celltype),
    !is.na(souporcell_origin)
  ) %>%
  mutate(
    classification = ifelse(
      celltype %in%
        c("HSC MPP", "MEP", "LMPP", "Cycling Progenitor", "Early GMP") &
        souporcell_origin == "recipient",
      yes = "recipient_progenitor",
      no = "other"
    )
  ) %>%
  group_by(patient_id, classification) %>%
  dplyr::count() %>%
  pivot_wider(names_from = classification, values_from = n) %>%
  mutate(percent = recipient_progenitor / (recipient_progenitor + other) * 100)
# --> "The persistent recipient HSPCs made up 0.27–1.31% of cells in these marrows"

# SOUPORCELL STATS ------------------------------------------------------------

# Figure 4A
seu@meta.data %>%
  tabyl("souporcell_origin") %>%
  adorn_totals(where = "row")
seu@meta.data %>%
  filter(!is.na(souporcell_origin)) %>%
  pull(patient_id) %>%
  unique %>%
  length
# Out of 33 patients in total, we successfully resolved cell origins in 21, yielding recipient/donor assignments for 342,780 out of 496,691 cells.

# NUMBAT CALLS ----------------------------------------------------------------

# Total percent normal and tumor
as_tibble(seu@meta.data) %>%
  pull(numbat_compartment) %>%
  tabyl()

# Select the same HSPC populations as in 3.3_DGE_HSPCs.R
as_tibble(seu@meta.data) %>%
  filter(
    !is.na(numbat_compartment),
    sample_status == "remission",
    celltype %in%
      c(
        "HSC MPP",
        "MEP",
        "LMPP",
        "Cycling Progenitor", # THIS IS A TYPO
        "Early GMP"
      )
  ) %>%
  dplyr::count(numbat_compartment, patient_id) %>%
  adorn_totals(where = "row")
# --> "In all patients where Numbat detected CNVs, every recipient HSPC in remission harbored them (75 cells across 3 patients)."

# TCR ENRICHMENT STATS --------------------------------------------------------

# Subset for T cells
seu_T <- subset(seu, !is.na(TCAT_Multinomial_Label))

# Some stats for the results section:
table(is.na(seu_T$CTstrict))
seu_T$CTstrict[!is.na(seu_T$CTstrict)] %>% unique %>% length()

# --> "We captured TCR sequences in 85.6% of all 199,957 T cells in our dataset, representing 115,205 unique clonotypes"
