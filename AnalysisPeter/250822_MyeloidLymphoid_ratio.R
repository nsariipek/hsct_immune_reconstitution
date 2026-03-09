# Peter van Galen, 250822
# Check myeloid/lymphoid ratio for Lev Silberstein

# Load libraries
library(Seurat)
library(tidyverse)
library(ggpubr)

# Set working directory
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/AnalysisPeter")

# Delete environment variables & load favorite function
rm(list = ls())

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Assign lineages
metadata <- as_tibble(seu@meta.data, rownames = "cell")

myeloid_celltypes <- c(
  "HSC MPP",
  "MEP",
  "EoBasoMast Precursor",
  "Megakaryocyte Precursor",
  "Early Erythroid",
  "Late Erythroid",
  "LMPP",
  "Cycling Progenitor",
  "Early GMP",
  "Late GMP",
  "Pro-Monocyte",
  "CD14 Mono",
  "CD16 Mono",
  "cDC",
  "pDC"
)

lymphoid_celltypes <- c(
  "Early Lymphoid",
  "Pro-B",
  "Pre-B",
  "Immature B",
  "Mature B",
  "Plasma Cell",
  "CD4 Naive",
  "CD4 Central Memory",
  "CD4 Effector Memory",
  "CD4 Regulatory",
  "CD8 Naive",
  "CD8 Central Memory",
  "CD8 Effector Memory 1",
  "CD8 Effector Memory 2",
  "CD8 Tissue Resident Memory",
  "T Proliferating",
  "NK",
  "NK CD56high",
  "NK Proliferating"
)

metadata <- metadata %>%
  mutate(
    myeloid_lymphoid = case_when(
      celltype %in% myeloid_celltypes ~ "myeloid",
      celltype %in% lymphoid_celltypes ~ "lymphoid",
      TRUE ~ NA
    )
  )

# Check
metadata %>% group_by(celltype, myeloid_lymphoid) %>% count %>% print(n = 40)

# Calculate proportion myeloid and lymphoid
metadata_summary <- metadata %>%
  filter(
    sample_status == "remission",
    timepoint %in% c(3, 5, 6),
    library_type == "MNC",
    #souporcell_origin == "donor",
    celltype != "Stromal"
  ) %>%
  drop_na(celltype) %>%
  group_by(patient_id, cohort, TP53_status, myeloid_lymphoid) %>%
  summarize(count = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(percent = count / sum(count) * 100)

# Stacked barplot
p1 <- metadata_summary %>%
  ggplot(aes(x = patient_id, y = percent, fill = myeloid_lymphoid)) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "Patient ID", y = "Percent of all cells") +
  theme_bw() +
  theme(
    aspect.ratio = 0.6,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black")
  )

# View & save
p1
ggsave("250822_1_Myeloid-Lymphoid_StackedBars.pdf", width = 8)


# Myeloid proportion by cohort
p2 <- metadata_summary %>%
  filter(myeloid_lymphoid == "myeloid") %>%
  ggplot(aes(x = cohort, y = percent)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = TP53_status), width = 0.2) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 90) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(y = "Myeloid cells (%)") +
  theme_bw() +
  theme(
    aspect.ratio = 2,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black")
  )

# View & save
p2
ggsave("250822_2_Myeloid_proportion_Cohort.pdf")

# Myeloid proportion by TP53 status
p3 <- metadata_summary %>%
  filter(myeloid_lymphoid == "myeloid") %>%
  mutate(TP53_status = factor(TP53_status, levels = c("WT", "MUT"))) %>%
  ggplot(aes(x = TP53_status, y = percent)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = patient_id), width = 0.2) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 90) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(y = "Myeloid cells (%)") +
  theme_bw() +
  theme(
    aspect.ratio = 2,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black")
  )

# View & save
p3
ggsave("250822_3_Myeloid_proportion_TP53.pdf")