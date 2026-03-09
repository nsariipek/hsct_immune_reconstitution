# Peter van Galen, 260110
# Plot candidate gene expression across datasets
# An updated version of this analysis is saved as 260110_Merge_and_annotate_three_studies.R and 260111_Expression_across_three_studies.R

library(tidyverse)
library(Seurat)
library(ggsci)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = T)
setwd(paste0(repo_root, "/AnalysisPeter/260127_R01_Figures"))

# Clear environment variables & load favorite function
rm(list = ls())
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Load data from three different projects—assuming the repositories are stored in the same folder
seu_aml_og <- readRDS("../../../reanalyze-aml2019/Seurat_AML.rds")
seu_bpdcn_og <- readRDS(
  "../../../Single-cell_BPDCN/04_XV-seq/BM_Seurat_Final.rds"
)
seu_hsct_og <- readRDS("../..//AuxiliaryFiles/250528_Seurat_complete.rds")

# Common genes
common_genes <- intersect(
  rownames(seu_aml_og),
  intersect(rownames(seu_bpdcn_og), rownames(seu_hsct_og))
)

# Get data from 2019 single-cell AML paper ------------------------------------

seu_aml <- subset(seu_aml_og, features = common_genes)
seu_aml <- CreateSeuratObject(
  counts = LayerData(seu_aml, assay = "RNA", layer = "counts"),
  meta.data = seu_aml@meta.data
)

# Filter for healthy control (BM) or malignant AML cells at Dx.
bm_ids <- seu_aml@meta.data |> filter(grepl("^BM", orig.ident)) |> rownames()
aml_ids <- seu_aml@meta.data |>
  filter(grepl("^AML.*D0", orig.ident), PredictionRefined == "malignant") |>
  rownames()
seu_aml <- subset(seu_aml, cells = c(bm_ids, aml_ids))

# Only keep necessary metadata columns
seu_aml@meta.data <- data.frame(
  orig.ident = seu_aml@meta.data$orig.ident,
  project = "aml_2019",
  tech = "Seq-well",
  donor_id = cutf(as.character(seu_aml@meta.data$orig.ident), d = "\\."),
  replicate = seu_aml@meta.data$orig.ident,
  cell_origin = case_when(
    grepl("BM", seu_aml@meta.data$orig.ident) ~ "Healthy",
    grepl("AML", seu_aml@meta.data$orig.ident) ~ "AML"
  ),
  select(seu_aml@meta.data, PredictionRefined),
  CellType1 = seu_aml@meta.data$CellType
) |>
  # Add leading 0 to BM IDs so they differ from the 2023 study
  mutate(donor_id = gsub("BM", "BM0", donor_id))

# Get healthy bone marrow donors from BPDCN project ---------------------------

seu_bpdcn <- subset(seu_bpdcn_og, features = common_genes)
seu_bpdcn <- CreateSeuratObject(
  counts = LayerData(seu_bpdcn, assay = "RNA", layer = "counts"),
  meta.data = seu_bpdcn@meta.data
)

# Only keep necessary metadata columns
seu_bpdcn@meta.data <- data.frame(
  orig.ident = seu_bpdcn@meta.data$orig.ident,
  project = "healthy_donors_2023",
  # fmt: skip
  tech = gsub("TenX", "10X_3prime",
    gsub("SW", "Seq-well", seu_bpdcn@meta.data$tech)
  ),
  donor_id = cutf(seu_bpdcn@meta.data$replicate, d = "\\."),
  replicate = seu_bpdcn@meta.data$replicate,
  cell_origin = "Healthy",
  CellType2 = seu_bpdcn@meta.data$CellType,
  row.names = colnames(seu_bpdcn)
) |>
  # Add five to donor ID so they differ from the 2019 study
  mutate(donor_id = sprintf("BM%02d", as.numeric(gsub("BM", "", donor_id)) + 5))

# Get data from 2026 HSCT paper -----------------------------------------------

seu_hsct <- subset(seu_hsct_og, features = common_genes)
seu_hsct <- CreateSeuratObject(
  counts = LayerData(seu_hsct, assay = "RNA", layer = "counts"),
  meta.data = seu_hsct@meta.data
)

# Filter for cells with celltype and souporcell assignments
exclude_ids <- seu_hsct_og@meta.data |>
  filter(is.na(celltype) | is.na(souporcell_origin)) |>
  rownames()
seu_hsct <- subset(seu_hsct, cells = setdiff(colnames(seu_hsct), exclude_ids))

# Also exclude early relapse cohort for better agreement with the HSCT paper
seu_hsct <- subset(seu_hsct, cohort_detail != "1-Early-relapse")

# Only keep necessary metadata columns
seu_hsct@meta.data <- data.frame(
  orig.ident = seu_hsct@meta.data$orig.ident,
  project = "hsct_2026",
  tech = "10X_5prime",
  donor_id = seu_hsct@meta.data$patient_id,
  select(
    seu_hsct@meta.data,
    sample_id,
    library_type,
    cohort,
    cohort_detail,
    timepoint,
    sample_status,
  ),
  cell_origin = case_when(
    seu_hsct@meta.data$souporcell_origin == "recipient" ~ "Recipient",
    seu_hsct@meta.data$souporcell_origin == "donor" ~ "Donor"
  ),
  numbat_compartment = seu_hsct@meta.data$numbat_compartment,
  CellType3 = seu_hsct@meta.data$celltype
)

# Merge and create consistent annotation --------------------------------------

# Merge three datasets together
seu_merge <- merge(seu_aml, merge(seu_bpdcn, seu_hsct))
seu_merge <- JoinLayers(seu_merge)
seu_merge <- NormalizeData(seu_merge)

# Add gene expression
metadata <- as_tibble(seu_merge@meta.data, rownames = "cell")
mygene <- "PRAME"
mygene <- "CALCRL"
mygene <- "RXFP1"
metadata$expr <- LayerData(seu_merge, layer = "data")[mygene, ]

# Create an overall cell type ID column. This is not the best approach: in future versions, I should apply BoneMarrowMap to aml_2019 and healthy_donors_2023, so that their annotations are more in line with hsct_2026
metadata <- metadata |>
  mutate(
    celltype = case_when(
      # aml_2019
      CellType1 == "HSC" ~ "HSPC",
      CellType1 == "Prog" ~ "HSPC",
      CellType1 == "earlyEry" ~ "earlyEry",
      CellType1 == "lateEry" ~ "lateEry",
      CellType1 == "GMP" ~ "GMP",
      CellType1 == "ProMono" ~ "ProMono",
      CellType1 == "Mono" ~ "Mono",
      CellType1 == "cDC" ~ "cDC",
      CellType1 == "pDC" ~ "pDC",
      CellType1 == "ProB" ~ "EarlyB",
      CellType1 == "Plasma" ~ "Plasma",
      CellType1 == "T" ~ "CD4",
      CellType1 == "CTL" ~ "CD8",
      CellType1 == "B" ~ "B",
      CellType1 == "NK" ~ "NK",
      CellType1 == "HSC-like" ~ "HSPC",
      CellType1 == "Prog-like" ~ "HSPC",
      CellType1 == "GMP-like" ~ "GMP",
      CellType1 == "ProMono-like" ~ "ProMono",
      CellType1 == "Mono-like" ~ "Mono",
      CellType1 == "cDC-like" ~ "cDC",
      # healthy_donors_2023
      CellType2 == "EarlyEry" ~ "EarlyEry",
      CellType2 == "Plasma" ~ "Plasma",
      CellType2 == "cDC" ~ "cDC",
      CellType2 == "HSPC" ~ "HSPC",
      CellType2 == "GMP" ~ "GMP",
      CellType2 == "LateEry" ~ "LateEry",
      CellType2 == "ProMono" ~ "ProMono",
      CellType2 == "ProB" ~ "EarlyB",
      CellType2 == "pDC" ~ "pDC",
      CellType2 == "Mono" ~ "Mono",
      CellType2 == "B" ~ "B",
      CellType2 == "CD8Memory" ~ "CD8",
      CellType2 == "NKT" ~ "CD8",
      CellType2 == "CD4Naive" ~ "CD4",
      CellType2 == "ncMono" ~ "Mono",
      CellType2 == "CD8Naive" ~ "CD8",
      CellType2 == "CD4Memory" ~ "CD4",
      CellType2 == "PreB" ~ "EarlyB",
      CellType2 == "CD8TermExh" ~ "CD8",
      CellType2 == "NK" ~ "NK",
      CellType2 == "GammaDeltaLike" ~ "CD8",
      # hsct_2026 - keep HSPC consistent with paper
      CellType3 == "HSC MPP" ~ "HSPC",
      CellType3 == "MEP" ~ "HSPC",
      CellType3 == "Cycling Progenitor" ~ "HSPC",
      CellType3 == "LMPP" ~ "HSPC",
      CellType3 == "Early GMP" ~ "HSPC", # Could also be GMP in this context
      CellType3 == "Late GMP" ~ "GMP",
      CellType3 == "Megakaryocyte Precursor" ~ "Prog", # Not HSPC in paper
      CellType3 == "EoBasoMast Precursor" ~ "Prog", # Not HSPC in paper
      CellType3 == "Early Lymphoid" ~ "Prog", # Not HSPC in paper
      CellType3 == "CD4 Effector Memory" ~ "CD4",
      CellType3 == "CD8 Tissue Resident Memory" ~ "CD8",
      CellType3 == "Plasma Cell" ~ "Plasma",
      CellType3 == "CD4 Naive" ~ "CD4",
      CellType3 == "CD4 Central Memory" ~ "CD4",
      CellType3 == "CD14 Mono" ~ "Mono",
      CellType3 == "NK" ~ "NK",
      CellType3 == "CD8 Naive" ~ "CD8",
      CellType3 == "Pro-Monocyte" ~ "ProMono",
      CellType3 == "CD8 Effector Memory 1" ~ "CD8",
      CellType3 == "cDC" ~ "cDC",
      CellType3 == "T Proliferating" ~ "CD8",
      CellType3 == "Early Erythroid" ~ "EarlyEry",
      CellType3 == "NK CD56high" ~ "NK",
      CellType3 == "CD8 Effector Memory 2" ~ "CD8",
      CellType3 == "Mature B" ~ "B",
      CellType3 == "Immature B" ~ "EarlyB",
      CellType3 == "CD8 Central Memory" ~ "CD8",
      CellType3 == "CD4 Regulatory" ~ "CD4",
      CellType3 == "Late Erythroid" ~ "LateEry",
      CellType3 == "pDC" ~ "pDC",
      CellType3 == "CD16 Mono" ~ "Mono",
      CellType3 == "Pre-B" ~ "EarlyB",
      CellType3 == "Pro-B" ~ "EarlyB",
      CellType3 == "NK Proliferating" ~ "NK",
      CellType3 == "Stromal" ~ "Stromal"
    )
  )

# Quick check
# select(metadata, CellType1, CellType2, CellType3, celltype) |>
#   unique |>
#   print(n = 80)

# Visualize (split by cell type) ----------------------------------------------

# Create tibble for plotting
plot_split_tib <- metadata |>
  filter(celltype %in% c("HSPC", "GMP", "ProMono", "Mono", "cDC")) |>
  # Remove non-malignant recipient cells from hsct_2026
  filter(is.na(numbat_compartment) | numbat_compartment != "normal") |>
  mutate(
    cell_origin = factor(
      cell_origin,
      levels = c("Healthy", "AML", "Donor", "Recipient")
    )
  ) |>
  mutate(CellOrigin_CellType = paste0(cell_origin, "_", celltype)) |>
  mutate(
    CellOrigin_CellType = factor(
      CellOrigin_CellType,
      levels = c(
        "Healthy_HSPC",
        "Healthy_GMP",
        "Healthy_ProMono",
        "Healthy_Mono",
        "Healthy_cDC",
        "AML_HSPC",
        "AML_GMP",
        "AML_ProMono",
        "AML_Mono",
        "AML_cDC",
        "Donor_HSPC",
        "Donor_GMP",
        "Donor_ProMono",
        "Donor_Mono",
        "Donor_cDC",
        "Recipient_HSPC",
        "Recipient_GMP",
        "Recipient_ProMono",
        "Recipient_Mono",
        "Recipient_cDC"
      )
    )
  ) |>
  mutate(
    project = case_when(
      project == "aml_2019" ~ "aml_diagnosis_2019",
      project == "healthy_donors_2023" ~ project,
      project == "hsct_2026" ~ paste0("hsct_", sample_status, "_2026"),
    )
  ) |>
  mutate(
    project = factor(
      project,
      levels = c(
        "healthy_donors_2023",
        "aml_diagnosis_2019",
        "hsct_pre-transplant_2026",
        "hsct_remission_2026",
        "hsct_relapse_2026"
      )
    )
  ) |>
  group_by(
    project,
    donor_id,
    cell_origin,
    celltype,
    CellOrigin_CellType
  ) |>
  filter(n() >= 10) |>
  summarize(
    mean_expr = mean(expr)
  )

# Wilcoxon tests. An alternative statistical approach would be DESeq2, but this is challenging if some datasets have a variable that is absent in others (such as donor/recipient)
wilcox_split_results <- plot_split_tib |>
  ungroup() |>
  filter(
    project %in%
      c("aml_diagnosis_2019", "hsct_remission_2026", "hsct_relapse_2026")
  ) |>
  summarize(
    # For the 2019 project, it will compare Healthy vs. AML. For the 2026 project, it will compare Donor vs. Recipient
    p_value = wilcox.test(
      x = mean_expr[cell_origin %in% c("Healthy", "Donor")],
      y = mean_expr[cell_origin %in% c("AML", "Recipient")]
    )$p.value,
    .by = c(celltype, project)
  ) |>
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    CellOrigin_CellType = ifelse(
      project == "aml_diagnosis_2019",
      paste0("AML_", celltype),
      paste0("Recipient_", celltype)
    )
  ) |>
  mutate(
    sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ "#",
      TRUE ~ ""
    )
  )

wilcox_split_results


# Visualize
plot_split_tib %>%
  ggplot(aes(x = CellOrigin_CellType, y = mean_expr)) +
  geom_boxplot(aes(fill = cell_origin), outlier.shape = NA, alpha = 0.3) +
  scale_fill_manual(
    values = c("grey", "firebrick", "darkgrey", "orange"),
    guide = guide_legend(ncol = 2)
  ) +
  geom_point(
    aes(color = donor_id),
    position = position_jitter(width = 0.3, height = 0)
  ) +
  scale_color_igv(guide = guide_legend(ncol = 4)) +
  geom_text(
    data = wilcox_split_results,
    aes(y = Inf, label = sig_label),
    vjust = 1.5,
    size = 4,
  ) +
  facet_grid(~project, scales = "free_x", space = "free_x") +
  labs(
    title = paste(mygene, "expression across three studies"),
    y = paste("Normalized", mygene, "expression")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(paste0("260110_", mygene, "_split_boxes.pdf"), width = 17, height = 5)

# Plot (not split by cell type) -----------------------------------------------

plot_tib <- metadata |>
  filter(celltype %in% c("HSPC", "GMP")) |>
  # Remove non-malignant recipient cells from hsct_2026
  filter(is.na(numbat_compartment) | numbat_compartment != "normal") |>
  mutate(
    cell_origin = factor(
      cell_origin,
      levels = c("Healthy", "AML", "Donor", "Recipient")
    )
  ) |>
  mutate(
    project = case_when(
      project == "aml_2019" ~ "aml_diagnosis_2019",
      project == "healthy_donors_2023" ~ project,
      project == "hsct_2026" ~ paste0("hsct_", sample_status, "_2026"),
    )
  ) |>
  mutate(
    project = factor(
      project,
      levels = c(
        "healthy_donors_2023",
        "aml_diagnosis_2019",
        "hsct_pre-transplant_2026",
        "hsct_remission_2026",
        "hsct_relapse_2026"
      )
    )
  ) |>
  group_by(
    project,
    donor_id,
    cell_origin
  ) |>
  filter(n() >= 10) |>
  summarize(
    mean_expr = mean(expr)
  )

# Wilcoxon tests
wilcox_results <- plot_tib |>
  ungroup() |>
  filter(
    project %in%
      c("aml_diagnosis_2019", "hsct_remission_2026", "hsct_relapse_2026")
  ) |>
  summarize(
    # For the 2019 project, it will compare Healthy vs. AML. For the 2026 project, it will compare Donor vs. Recipient
    p_value = wilcox.test(
      x = mean_expr[cell_origin %in% c("Healthy", "Donor")],
      y = mean_expr[cell_origin %in% c("AML", "Recipient")]
    )$p.value,
    .by = c(project)
  ) |>
  mutate(cell_origin = c("AML", "Recipient", "Recipient"))

wilcox_results


# Visualize
plot_tib %>%
  ggplot(aes(x = cell_origin, y = mean_expr)) +
  geom_boxplot(aes(fill = cell_origin), outlier.shape = NA, alpha = 0.3) +
  scale_fill_manual(
    values = c("grey", "firebrick", "darkgrey", "orange"),
    guide = guide_legend(ncol = 2)
  ) +
  geom_point(
    aes(color = donor_id),
    position = position_jitter(width = 0.3, height = 0)
  ) +
  scale_color_igv(guide = guide_legend(ncol = 4)) +
  geom_text(
    data = wilcox_results,
    aes(
      y = Inf,
      label = ifelse(
        p_value < 0.1,
        yes = format(p_value, scientific = TRUE, digits = 3),
        no = "n.s."
      )
    ),
    vjust = 1.5,
    size = 3,
  ) +
  #facet_grid(~project, scales = "free_x", space = "free_x") +
  facet_grid(
    ~project,
    scales = "free_x",
    space = "free_x",
    labeller = labeller(project = ~ gsub("_", "\n", .x))
  ) +
  labs(
    title = paste(mygene, "expression across three studies"),
    y = paste("Normalized", mygene, "expression")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11)
  )

ggsave(paste0("260110_", mygene, "_boxes.pdf"), width = 9, height = 5)
