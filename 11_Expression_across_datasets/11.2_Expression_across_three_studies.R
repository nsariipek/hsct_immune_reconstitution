# Peter van Galen, 260111
# Visualize gene expression across three studies

library(tidyverse)
library(Seurat)
library(janitor)
library(ggsci)
library(ggnewscale)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/AnalysisPeter"))

# Clear environment variables & load favorite function
rm(list = ls())
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Load Seurat data from three studies. This is saved in Terra's Google bucket and available upon request (see 11.1_Merge_and_annotate_three_studies.R).
seu <- readRDS("../AuxiliaryFiles/260110_Seurat_ThreeStudies.rds")


# Plot ------------------------------------------------------------------------

# Add gene expression
metadata <- as_tibble(seu@meta.data, rownames = "cell")
mygene <- "PRAME"
mygene <- "CALCRL"
mygene <- "RXFP1"
metadata$expr <- LayerData(seu, layer = "data")[mygene, ]

# Create tibble for plotting
wrangle_tib <- metadata |>
  # Consolidate some cell types for clarity of box plots
  mutate(
    celltype_tmp = case_when(
      celltype %in%
        c("HSC MPP", "MEP", "LMPP", "Cycling Progenitor", "Early GMP") ~ "HSPC",
      celltype == "Late GMP" ~ "LateGMP",
      celltype %in% c("Pro-Monocyte", "CD14 Mono", "CD16 Mono") ~ "Mono",
      .default = celltype
    )
  ) |>
  # Filter for cell types most likely to contain malignant AML cells (`metadata |> filter(PredictionRefined == "malignant") |> pull(celltype) |> tabyl() |> mutate(percent = percent*100) |> filter(n > 100)`)
  filter(celltype_tmp %in% c("HSPC", "LateGMP", "Mono")) |>
  # Remove non-malignant recipient cells from hsct_2026
  filter(is.na(numbat_compartment) | numbat_compartment != "normal") |>
  mutate(
    cell_origin = factor(
      cell_origin,
      levels = c("Healthy", "AML", "Donor", "Recipient")
    )
  ) |>
  mutate(CellOrigin_CellType = paste0(cell_origin, "_", celltype_tmp)) |>
  mutate(
    CellOrigin_CellType = factor(
      CellOrigin_CellType,
      levels = c(
        "Healthy_HSPC",
        "Healthy_LateGMP",
        "Healthy_Mono",
        "AML_HSPC",
        "AML_LateGMP",
        "AML_Mono",
        "Donor_HSPC",
        "Donor_LateGMP",
        "Donor_Mono",
        "Recipient_HSPC",
        "Recipient_LateGMP",
        "Recipient_Mono"
      )
    )
  ) |>
  mutate(
    project = case_when(
      project == "aml_2019" ~ "aml_2019_diagnosis",
      project == "healthy_donors_2023" ~ project,
      project == "hsct_2026" ~ paste0(project, "_", sample_status),
    )
  ) |>
  mutate(
    project = factor(
      project,
      levels = c(
        "aml_2019_diagnosis",
        "healthy_donors_2023",
        "hsct_2026_pre-transplant",
        "hsct_2026_remission",
        "hsct_2026_relapse"
      )
    )
  )

# Summarize to get mean expression
plot_tib <- wrangle_tib |>
  group_by(
    project,
    donor_id,
    cell_origin,
    celltype_tmp,
    CellOrigin_CellType
  ) |>
  filter(n() >= 10) |>
  summarize(mean_expr = mean(expr))

# Wilcoxon tests. An alternative statistical approach would be DESeq2, but this is challenging if some datasets have a variable that is absent in others (such as donor/recipient)
wilcox_results <- plot_tib |>
  ungroup() |>
  filter(
    project %in%
      c("aml_2019_diagnosis", "hsct_2026_remission", "hsct_2026_relapse")
  ) |>
  summarize(
    # For the 2019 project, it will compare Healthy vs. AML. For the 2026 project, it will compare Donor vs. Recipient
    p_value = wilcox.test(
      mean_expr[cell_origin %in% c("Healthy", "Donor")],
      mean_expr[cell_origin %in% c("AML", "Recipient")]
    )$p.value,
    .by = c(celltype_tmp, project)
  ) |>
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    CellOrigin_CellType = ifelse(
      project == "aml_2019_diagnosis",
      paste0("AML_", celltype_tmp),
      paste0("Recipient_", celltype_tmp)
    )
  ) |>
  mutate(
    sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

wilcox_results

# Visualize
plot_tib %>%
  ggplot(aes(x = CellOrigin_CellType, y = mean_expr)) +
  geom_boxplot(aes(fill = cell_origin), outlier.shape = NA, alpha = 0.3) +
  scale_fill_manual(
    values = c("grey", "firebrick", "darkgrey", "orange"),
    guide = guide_legend(ncol = 2)
  ) +
  new_scale_fill() +
  geom_point(
    aes(fill = donor_id),
    shape = 21,
    size = 2,
    color = "black",
    stroke = 0.3,
    position = position_jitter(width = 0.3, height = 0)
  ) +
  scale_fill_igv(guide = guide_legend(ncol = 4)) +
  geom_text(
    data = wilcox_results,
    aes(y = Inf, label = sig_label),
    vjust = 1.5,
    size = 4,
  ) +
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
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(paste0("11.2_", mygene, "_split_boxes.pdf"), width = 12, height = 5)


# Plot (not split by cell type) -----------------------------------------------

plot2_tib <- wrangle_tib |>
  group_by(
    project,
    donor_id,
    cell_origin
  ) |>
  filter(n() >= 10) |>
  summarize(mean_expr = mean(expr))

# Wilcoxon tests
wilcox_results2 <- plot2_tib |>
  ungroup() |>
  filter(
    project %in%
      c("aml_2019_diagnosis", "hsct_2026_remission", "hsct_2026_relapse")
  ) |>
  summarize(
    # For the 2019 project, it will compare Healthy vs. AML. For the 2026 project, it will compare Donor vs. Recipient
    p_value = wilcox.test(
      x = mean_expr[cell_origin %in% c("Healthy", "Donor")],
      y = mean_expr[cell_origin %in% c("AML", "Recipient")]
    )$p.value,
    .by = c(project)
  )

wilcox_results2


# Visualize
plot2_tib %>%
  ggplot(aes(x = cell_origin, y = mean_expr)) +
  geom_boxplot(aes(fill = cell_origin), outlier.shape = NA, alpha = 0.3) +
  scale_fill_manual(
    values = c("grey", "firebrick", "darkgrey", "orange"),
    guide = guide_legend(ncol = 2)
  ) +
  new_scale_fill() +
  geom_point(
    aes(fill = donor_id),
    shape = 21,
    size = 2,
    color = "black",
    stroke = 0.3,
    position = position_jitter(width = 0.3, height = 0)
  ) +
  scale_fill_igv(guide = guide_legend(ncol = 4)) +
  geom_text(
    data = wilcox_results2,
    aes(
      x = 1.5,
      y = Inf,
      label = ifelse(
        p_value < 0.05,
        yes = paste0("P=", format(p_value, scientific = TRUE, digits = 3)),
        no = "n.s."
      )
    ),
    inherit.aes = FALSE,
    vjust = 1.5,
    size = 3,
  ) +
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

ggsave(paste0("11.2_", mygene, "_boxes.pdf"), width = 9, height = 5)
