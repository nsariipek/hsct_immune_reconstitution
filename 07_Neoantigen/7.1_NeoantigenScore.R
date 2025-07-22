# Nurefsan Sariipek and Peter van Galen, 250719
# Quantify neoantigen score per clonotype

# Load the libraries
library(tidyverse)
library(Seurat)
library(ggforce)
library(ggnewscale)
library(ggrastr)

# Start with a clean slate
rm(list = ls())

# Set working directory and load Seurat object
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/07_Neoantigen")
# For VM:
#setwd("~/TP53_ImmuneEscape/10_Neoantigen/")
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Subset for T cells
seu_T <- subset(seu, subset = !is.na(TCAT_Multinomial_Label))

# Add ASA Score from starCAT (Kotliar 2025)
usage_tib <- read_tsv("../05_DGE/5.1_starCAT/starCAT.scores.txt.gz") %>%
  column_to_rownames("...1")
seu_T <- AddMetaData(seu_T, select(usage_tib, ASA))


#### Score signatures ####

# Load antigen recognition signatures from Rosenberg lab: Lowery 2022, Table S10 (https://www.sciencedirect.com/science/article/pii/S1535610823003963)
cd4neoA <- read.csv("signatures/cd4.csv")
cd8neoA <- read.csv("signatures/cd8.csv")
# And Yossef 2023 (https://www.science.org/doi/10.1126/science.abl5447)
neoA <- read.csv("signatures/neoantigen.csv")

# Subset for CD8 T cell subsets with (similar to 4.1_Monocle3.R). Excluding naive cells changes the results but overall conclusions stay the same.
seu_subset <- subset(
  seu_T,
  subset = TCAT_Multinomial_Label %in%
    c("CD8_Naive", "CD8_CM", "CD8_EM", "CD8_TEMRA") &
    timepoint %in% c(3, 5, 6) &
    sample_status == "remission"
)
T_marker <- "CD8"

# Alternatively, subset for CD4 T cell subsets (similar to 4.1).
seu_subset <- subset(
  seu_T,
  subset = TCAT_Multinomial_Label %in%
    c("CD4_Naive", "CD4_CM", "CD4_EM") &
    timepoint %in% c(3, 5, 6) &
    sample_status == "remission"
)
T_marker <- "CD4"

# Add module scores to the subsetted dataset
seu_subset <- NormalizeData(seu_subset, assay = "RNA")
seu_subset <- AddModuleScore(
  seu_subset,
  features = list(cd8neoA$gene, neoA$gen, cd4neoA$gene),
  name = c("cd4neoA", "cd8neoA", "neoA"),
  search = T
)


#### WRANGLE ####

# Select one. Usage is based on the original papers, cited above.
current_sig <- "ASA" # can be used for both
current_sig <- "cd4neoA1" # only use for CD4 T subset
current_sig <- "cd8neoA2" # only use for CD8 subset
current_sig <- "neoA3" # only use for CD8 subset

# Get metadata & select needed variables
metadata <- as_tibble(seu_subset@meta.data) %>%
  select(
    patient_id,
    CTstrict,
    cohort,
    patient_id,
    TP53_status,
    all_of(current_sig)
  ) %>%
  mutate(TP53_status = factor(TP53_status, levels = c("WT", "MUT")))

# For each clonotype, calculate relative size and median signature score. We don't actually plot the size (but used to in a prior version)
ct_tib <- metadata %>%
  group_by(patient_id) %>%
  mutate(n_total = n(), .groups = "drop") %>%
  group_by(patient_id, cohort, TP53_status, CTstrict) %>%
  summarize(
    n = n(),
    prop = n / unique(n_total),
    ct_median = median(.data[[current_sig]]),
  )

# For each patient, calculate median signature score
pt_tib <- ct_tib %>%
  group_by(patient_id, TP53_status) %>%
  summarize(pt_median = median(ct_median))


#### STATS ####

# Statistics v1: use clonotypes as replicates (too lenient)
ct_MUT_meanScores <- filter(ct_tib, TP53_status == "MUT") %>% pull(ct_median)
ct_WT_meanScores <- filter(ct_tib, TP53_status == "WT") %>% pull(ct_median)
wilcox.test(ct_MUT_meanScores, ct_WT_meanScores)

# Statistics v2: use patients as replicates (ok)
pt_MUT_meanScores <- filter(pt_tib, TP53_status == "MUT") %>% pull(pt_median)
pt_WT_meanScores <- filter(pt_tib, TP53_status == "WT") %>% pull(pt_median)
wcox <- wilcox.test(pt_MUT_meanScores, pt_WT_meanScores)
# Fold change
fc <- median(pt_MUT_meanScores) / median(pt_WT_meanScores)


#### PLOT ####

# Order violins
pt_high2low <- as.character(
  pt_tib %>% arrange(desc(pt_median)) %>% pull(patient_id)
)
ct_tib <- ct_tib %>%
  mutate(patient_id = fct_relevel(patient_id, pt_high2low))

ct_tib %>%
  ggplot(aes(x = patient_id, y = ct_median)) +
  geom_tile(
    data = distinct(select(ct_tib, patient_id, cohort, TP53_status)),
    aes(y = max(ct_tib$ct_median) + 0.06, fill = cohort),
    height = 0.03,
    width = 0.9
  ) +
  scale_fill_manual(
    values = c(`long-term-remission` = "#4775FF", relapse = "#E64B35")
  ) +
  new_scale_fill() +
  geom_tile(
    data = distinct(select(ct_tib, patient_id, cohort, TP53_status)),
    aes(y = max(ct_tib$ct_median) + 0.02, fill = TP53_status),
    height = 0.03,
    width = 0.9,
    show.legend = F
  ) +
  rasterize(
    geom_sina(aes(color = TP53_status), shape = 16, size = 1, scale = "width"),
    dpi = 600
  ) +
  geom_violin(alpha = 0, scale = "width", draw_quantiles = 0.5) +
  scale_color_manual(values = c(WT = "#7FC97F80", MUT = "#66009980")) +
  scale_fill_manual(values = c(WT = "#7FC97F", MUT = "#660099")) +
  theme_bw() +
  theme(
    aspect.ratio = 0.5,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save as pdf
ggsave(
  paste0("7.1_Ordered_violins_", T_marker, "_", current_sig, ".pdf"),
  width = 8,
  height = 4
)

pt_tib %>%
  ggplot(aes(x = TP53_status, y = pt_median, fill = TP53_status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = patient_id), width = 0.2, size = 2) +
  scale_fill_manual(values = c(WT = "#7FC97F", MUT = "#660099")) +
  annotate(
    "text",
    x = 1.5,
    y = max(pt_tib$pt_median) - 0.01,
    label = paste0(
      "Fold change: ",
      round(fc, 4),
      "\n",
      "p = ",
      signif(wcox$p.value, digits = 4)
    )
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 2,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  )

# Save as pdf
ggsave(
  paste0("7.2_Patient_boxplot_", T_marker, "_", current_sig, ".pdf"),
  width = 4,
  height = 5
)
