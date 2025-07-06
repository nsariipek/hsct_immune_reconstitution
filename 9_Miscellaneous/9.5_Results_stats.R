# Peter van Galen, 250705
# Extract metrics for the paper results section

# Load libraries
library(Seurat)
library(tidyverse)
library(janitor)

# Set working directory
setwd(
  "~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Miscellaneous"
)

# Delete environment variables & load favorite function
rm(list = ls())
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")


### TCR ENRICHMENT ###

# Subset for T cells
T_celltypes <- c(
  "CD4 Naive",
  "CD4 Central Memory",
  "CD4 Effector Memory",
  "CD4 Regulatory",
  "CD8 Naive",
  "CD8 Central Memory",
  "CD8 Effector Memory 1",
  "CD8 Effector Memory 2",
  "CD8 Tissue Resident Memory",
  "T Proliferating"
)
seu_T <- subset(seu, subset = celltype %in% T_celltypes)

# Some stats for the results section:
table(is.na(seu_T$CTstrict))
seu_T$CTstrict[!is.na(seu_T$CTstrict)] %>% unique %>% length()

# --> "We captured TCR sequences in 85.6% of all 199,957 T cells in our dataset, representing 115,205 unique clonotypes"

### SOUPORCELL CALLS ###

# Number of donor and recipient HSPCs at remission
as_tibble(seu@meta.data) %>%
  filter(
    sample_status == "remission",
    celltype %in%
      c(
        "HSC MPP",
        "MEP",
        "LMPP",
        "Cycling Progenitors",
        "Early GMP"#,
        #"Late GMP"
      )
  ) %>%
  pull(souporcell_origin) %>%
  tabyl()
# --> "Our dataset provides a unique opportunity to compare recipient HSPCs (n=525 cells) to donor HSPCs (n=1,466 cells) at remission"

# Similar to 3_DGE/3.3_DGE_progenitors.R
# TODO: CHECK BY RERUNNING THAT SCRIPT!
as_tibble(seu@meta.data) %>%
  filter(
cohort == "relapse",
timepoint %in% c("3", "5", "6"),
sample_status == "remission",
    celltype %in%
      c(
        "HSC MPP",
        "MEP",
        "LMPP",
        "Cycling Progenitors",
        "Early GMP"
      ),
      !is.na(souporcell_origin)
  ) %>%
  pull(patient_id) %>% as.character %>% unique %>% sort
# --> "For this analysis, we excluded the long-term remission cohort who did not have persistent recipient HSPCs and P30-P33 who did not have 3-6 month remission samples, leaving seven patients from the relapse cohort."

### NUMBAT CALLS ###

# Total percent normal and tumor
as_tibble(seu@meta.data) %>%
  pull(numbat_compartment) %>%
  tabyl()

# Total percent normal and tumor
as_tibble(seu@meta.data) %>%
  filter(sample_status == "remission") %>%
  pull(numbat_compartment) %>%
  tabyl()
# --> "At clinical remission time points, 14.0% of recipient cells were classified as malignant"

# Select the same HSPC populations as in 6.4_Plots.R
as_tibble(seu@meta.data) %>%
  filter(
    sample_status == "remission",
    !is.na(numbat_compartment),
    celltype %in%
      c(
        "HSC MPP",
        "MEP",
        "LMPP",
        "Cycling Progenitors",
        "Early GMP"
      )
  ) %>%
  group_by(numbat_compartment, patient_id) %>%
  dplyr::count()
# --> "Strikingly, 100% of recipient HSPCs were classified as malignant (58 cells across 3 patients)."
