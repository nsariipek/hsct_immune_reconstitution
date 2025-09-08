# Nurefsan Sariipek and Peter van Galen, 250718
# Analyze diversity in remission samples ~3 months after transplant using subsampling based on cell numbers which is different from scRepertoire built-in subsampling

# Load the libraries
library(tidyverse)
library(Seurat)
library(ggpubr)

# Empty environment
rm(list = ls())

# Set working directory for Nurefsan/Terra:
setwd("~/hsct_immune_reconstitution/06_TCR_Diversity/")

# For Peter:
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/06_TCR_Diversity")

# Load final Seurat object including TCR calls & subset for T cells
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")
seu_T <- subset(seu, subset = !is.na(TCAT_Multinomial_Label))

# Load colors
celltype_colors_df <- read.table(
  "../celltype_colors.txt",
  sep = "\t",
  header = T,
  stringsAsFactors = F,
  comment.char = ""
)
celltype_colors <- setNames(
  celltype_colors_df$color,
  celltype_colors_df$celltype
)
cohort_colors <- c("long-term-remission" = "#546fb5", "relapse" = "#e54c35")


######## PREPARE DATA ########

# Turn metadata to a tibble and keep only needed variables
metadata_tib <- as_tibble(seu_T@meta.data, rownames = "cell") %>%
  select(
    cell,
    cohort,
    patient_id,
    timepoint,
    sample_status,
    TP53_status,
    TCAT_Multinomial_Label,
    CTstrict,
    souporcell_origin
  )

# Filter for remission samples ~3M after transplant
metasubset1_tib <- metadata_tib %>%
  filter(
    sample_status == "remission",
    timepoint %in% c(3, 5, 6)
  )


######## SUBSET ########

# Subsetting for different cell types. Choose one!

# All cells
cell_subset <- "_all"
metasubset2_tib <- metasubset1_tib

# CD4+ T cells
cell_subset <- "_CD4"
metasubset2_tib <- metasubset1_tib %>%
  subset(TCAT_Multinomial_Label %in% c("CD4_Naive", "CD4_CM", "CD4_EM", "Treg"))

# CD8+ T cells
cell_subset <- "_CD8"
metasubset2_tib <- metasubset1_tib %>%
  subset(
    TCAT_Multinomial_Label %in% c("CD8_Naive", "CD8_CM", "CD8_EM", "CD8_TEMRA")
  )

# Recipient cells
cell_subset <- "_recipient"
metasubset2_tib <- metasubset1_tib %>% subset(souporcell_origin == "recipient")

# Donor cells
cell_subset <- "_donor"
metasubset2_tib <- metasubset1_tib %>% subset(souporcell_origin == "donor")


######## CALCULATE AND PLOT DIVERSITY ########

# Cell threshold
min_cells <- 300

# Visualize number of cells per patient
metasubset2_tib %>%
  group_by(patient_id) %>%
  count() %>%
  ggplot(aes(x = patient_id, y = n)) +
  geom_point() +
  geom_hline(yintercept = min_cells, color = "red") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# Filter for samples with sufficient cells
metasubset3_tib <- metasubset2_tib %>%
  group_by(patient_id) %>%
  filter(n() >= min_cells) %>%
  mutate(patient_id = as.character(patient_id))

# Sample min_cells from each patient. This gives all the plots below the same y-axis range
set.seed(94)
metasubset3_tib <- metasubset3_tib %>%
  group_by(patient_id) %>%
  slice_sample(n = min_cells)

# Load diversity functions, split data into a list, and calculate diversity
source("DiversityFunctions.R")
metasubset_ls <- split(metasubset3_tib, f = metasubset3_tib$patient_id)
diversities_df <- compute_diversity(metasubset_ls, "CTstrict", 1000)
diversities_df <- diversities_df %>% rename(patient_id = "group")

# Add information to make annotated plots
metasubset_summary <- metasubset3_tib %>%
  group_by(patient_id, cohort, TP53_status) %>%
  summarize
joined_tibble <- diversities_df %>% left_join(metasubset_summary)

# Visualize barplot
plot_TCR <- function(df, title = NULL) {
  df %>%
    ggplot(aes(x = cohort, y = inv.simpson, fill = cohort)) +
    geom_boxplot(width = 0.7, alpha = 0.9, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 2, color = "black") +
    scale_fill_manual(values = cohort_colors) +
    coord_cartesian(ylim = c(0, 310)) +
    labs(y = "Inverse Simpson Index") +
    stat_compare_means(
      aes(group = cohort),
      method = "wilcox.test",
      label = "p.format",
      label.y = 300
    ) +
    theme_minimal() +
    theme(
      aspect.ratio = 2,
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(
        color = "black",
        angle = 45,
        hjust = 1
      ),
      axis.title.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.ticks = element_line(color = "black"),
      legend.position = "none"
    )
}

TCR_all <- plot_TCR(joined_tibble, title = "All samples")
TCR_MT <- plot_TCR(
  filter(joined_tibble, TP53_status == "MUT"),
  title = "TP53-MUT samples"
)
TCR_WT <- plot_TCR(
  filter(joined_tibble, TP53_status == "WT"),
  title = "TP53-WT samples"
)

# Save as a pdf: combined TP53-MUT and TP53-WT, or separate
pdf(
  paste0("6.4_TCR_Diversity_3M", cell_subset, "_both.pdf"),
  width = 2.5,
  height = 3.5
)
TCR_all
dev.off()

pdf(
  paste0("6.4_TCR_Diversity_3M", cell_subset, "_MT.pdf"),
  width = 2.5,
  height = 3.5
)
TCR_MT
dev.off()

pdf(
  paste0("6.4_TCR_Diversity_3M", cell_subset, "_WT.pdf"),
  width = 2.5,
  height = 3.5
)
TCR_WT
dev.off()


######## MISCELLANEOUS THOUGHT ########

# It may appear surprising that recipient T cells appear more diverse than donor T cells,
# but it makes sense considering a much higher fraction of recipient T cells are CD4+:
metasubset1_tib %>%
  filter(!is.na(souporcell_origin)) %>%
  dplyr::count(souporcell_origin, TCAT_Multinomial_Label) %>%
  pivot_wider(names_from = souporcell_origin, values_from = n) %>%
  mutate(
    fraction_recipient = recipient / sum(recipient),
    fraction_donor = donor / sum(donor)
  ) %>%
  pivot_longer(cols = c(fraction_recipient, fraction_donor)) %>%
  ggplot(aes(x = name, y = value, fill = TCAT_Multinomial_Label)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = celltype_colors) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )
