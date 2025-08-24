# Peter van Galen, 250705
# Extract metrics for the paper results section

# Load libraries
library(Seurat)
library(tidyverse)
library(janitor)

# Set working directory
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/10_Miscellaneous")

# Delete environment variables & load favorite function
rm(list = ls())
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")


### TCR ENRICHMENT ###

# Subset for T cells
seu_T <- subset(seu, !is.na(TCAT_Multinomial_Label))

# Some stats for the results section:
table(is.na(seu_T$CTstrict))
seu_T$CTstrict[!is.na(seu_T$CTstrict)] %>% unique %>% length()

# --> "We captured TCR sequences in 85.6% of all 199,957 T cells in our dataset, representing 115,205 unique clonotypes"

### SOUPORCELL CALLS ###

# Figure 4A
seu@meta.data %>% tabyl("souporcell_origin") %>% adorn_totals(where = "row")
seu@meta.data %>%
  filter(!is.na(souporcell_origin)) %>%
  pull(patient_id) %>%
  unique %>%
  length
seu@meta.data %>% tabyl("numbat_compartment") %>% adorn_totals(where = "row")
seu@meta.data %>%
  filter(!is.na(numbat_compartment)) %>%
  pull(patient_id) %>%
  unique %>%
  length

# Recipient and donor HSPCs at remission. See 05_DGE/5.3_DGE_HSPCs.R for more details
meta_subset <- as_tibble(seu@meta.data) %>%
  filter(
    sample_status == "remission",
    timepoint %in% c(3, 5, 6),
    celltype %in%
      c("HSC MPP", "MEP", "LMPP", "Cycling Progenitors", "Early GMP"),
    souporcell_origin %in% c("donor", "recipient"),
    cohort == "relapse",
    sample_id != "P23_Rem1" # does not have recipient HSPCs
  )
meta_subset %>%
  pull(patient_id) %>%
  as.character %>%
  unique %>%
  sort
# --> "For this analysis, we excluded the long-term remission cohort which had <10 persistent recipient HSPCs and P30-P33 who did not have 3-6 month remission samples, leaving six patients from the relapse cohort."

meta_subset %>%
  pull(souporcell_origin) %>%
  tabyl()
# --> "Genes that were upregulated in recipient HSPCs (n=417 cells) compared to their donor counterparts (n=439) included"

# Similar to 8.4_Plots.R
merged_prog_counts_tib <- as_tibble(seu@meta.data) %>%
  filter(
    sample_status == "remission",
    timepoint %in% c(3, 5, 6),
    celltype %in%
      c("HSC MPP", "MEP", "LMPP", "Cycling Progenitors", "Early GMP")
  ) %>%
  dplyr::count(patient_id, cohort, souporcell_origin)
# Pivot wider to compute proportion donor
merged_prog_counts_tib %>%
  filter(souporcell_origin %in% c("donor", "recipient")) %>%
  pivot_wider(
    names_from = souporcell_origin,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(donor_percentage = donor / (donor + recipient) * 100) %>%
  arrange(donor_percentage) %>%
  print(n = 30)
# --> "All patients in whom recipient HSPCs persisted at >5% developed relapse within 18 months"

### NUMBAT CALLS ###

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
        "Cycling Progenitors",
        "Early GMP"
      )
  ) %>%
  dplyr::count(numbat_compartment, patient_id) %>%
  adorn_totals(where = "row")
# --> "In all patients where Numbat detected CNVs, every recipient HSPC in remission harbored them (58 cells across 3 patients)."
