# Peter van Galen, 250418
# Quickly check T cell type proportions between cohorts for the TP53 mutated patients

# Clear environment
rm(list = ls())

# Load complete Seurat object
seu <- readRDS("../AuxiliaryFiles/250418_Seurat_all_cells_annotated.rds")

# Define T and NK cell types
t_celltypes <- c("CD4 Naive", "CD4 Memory", "CD4 Effector Memory", "Treg",
"CD8 Naive", "CD8 Memory", "CD8 Effector", "CD8 Exhausted",
"Gamma-Delta T")

# Quick check of TP53 cell type proportions
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell") 
celltype_proportions <- metadata_tib %>%
  filter(TP53_status == "MUT", timepoint %in% c(3, 5, 6),
        sample_status == "remission", celltype %in% t_celltypes) %>%
  group_by(cohort, patient_id, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cohort, patient_id) %>%
  mutate(percent = count / sum(count) * 100)

celltype_proportions %>%
  ggplot(aes(x = cohort, y = percent)) +
  geom_boxplot(outliers = F) +
  geom_jitter(aes(color = patient_id)) +
  facet_wrap(~celltype, ncol = 4) + #, scales = "free_y"
  stat_compare_means(method = "wilcox.test",
                             label = "p.format",
                             vjust = 1.5) +
  theme_bw() +
  theme(panel.grid = element_blank())


# Try with TCAT annotations
usage_tib <- read_tsv("../3_DGE/3.1_starCAT/starCAT.scores.txt.gz") %>% rename("cell" = "...1")
tcat_celltypes_df <- data.frame(select(usage_tib, Multinomial_Label), row.names = usage_tib$cell)
seu <- AddMetaData(seu, tcat_celltypes_df)
seu$Multinomial_Label <- factor(seu$Multinomial_Label, levels = c("CD4_Naive", "CD4_CM", "CD4_EM", "Treg", "CD8_Naive", "CD8_CM", "CD8_EM",  "CD8_TEMRA", "MAIT", "gdT"))

# Quick check of TP53 cell type proportions
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell") 
celltype_proportions <- metadata_tib %>%
  filter(TP53_status == "MUT", timepoint %in% c(3, 5, 6),
        sample_status == "remission", celltype %in% t_celltypes) %>%
  group_by(cohort, patient_id, Multinomial_Label) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cohort, patient_id) %>%
  mutate(percent = count / sum(count) * 100)

celltype_proportions %>%
  ggplot(aes(x = cohort, y = percent)) +
  geom_boxplot(outliers = F) +
  geom_jitter(aes(color = patient_id)) +
  facet_wrap(~Multinomial_Label) + #, scales = "free_y"
  stat_compare_means(method = "wilcox.test",
                             label = "p.format",
                             vjust = 1.5) +
  theme_bw() +
  theme(panel.grid = element_blank())
