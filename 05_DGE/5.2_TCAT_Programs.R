# Peter van Galen and Nurefsan Sariipek, 250716
# Assess TCAT programs in T cells

library(Seurat)
library(tidyverse)
library(scattermore)
library(ggpubr)
library(ggforce)
library(ComplexHeatmap)
library(circlize)

# Set working directory
#VM: setwd("~/hsct_immune_reconstitution/05_DGE/")
#Local:
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/05_DGE")

# Delete environment variables & load favorite function
rm(list = ls())
cutf <- function(x, f = 1, d = "/") {
      sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")
seu_T <- subset(seu, !is.na(TCAT_Multinomial_Label))

# Load colors from 2.3_PvG-Colors.R
celltype_colors_df <- read.table(
      "../celltype_colors.txt",
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      comment.char = ""
)
celltype_colors <- setNames(
      celltype_colors_df$color,
      celltype_colors_df$celltype
)

# Combine TCAT scores and program usage results with Seurat metadata
usage_tib <- read_tsv("5.1_starCAT/starCAT.scores.txt.gz") %>%
      rename("cell" = "...1")
scores_tib <- read_tsv("5.1_starCAT/starCAT.rf_usage_normalized.txt.gz") %>%
      rename("cell" = "...1")
metadata_tib <- as_tibble(seu_T@meta.data, rownames = "cell")
metadata_tib <- left_join(metadata_tib, scores_tib)
metadata_tib <- left_join(metadata_tib, usage_tib)

# Filter data for remission samples ~3M after transplant
metadata_subset <- metadata_tib %>%
      filter(
            sample_status == "remission",
            timepoint %in% c(3, 5, 6)
      )

# Subset for analysis of activity programs from the bioRxiv *CAT paper (https://doi.org/10.1101/2024.05.03.592310, Figure 2C, Table S1)
activity_programs <- c(
      "CellCycle-S",
      "CellCycle-Late-S",
      "CellCycle-G2M",
      "BCL2/FAM13A",
      "ICOS/CD38",
      "CD172a/MERTK",
      "SOX4/TOX2",
      "Heatshock",
      "ISG",
      "OX40/EBI3",
      "CD40LG/TXNIP",
      "Cytoskeleton",
      "Cytotoxic",
      "CTLA4/CD38",
      "Exhaustion",
      "HLA",
      "NME1/FABP5",
      "Translation",
      "IEG",
      "TIMD4/TIM3",
      "IL10/IL19",
      "RGCC/MYADM",
      "Metallothionein",
      "Multi-Cytokine"
)
identity_programs <- c(
      "CD8-Naive",
      "CD8-EM",
      "CD8-Trm",
      "TEMRA",
      "CD4-Naive",
      "CD4-CM",
      "Treg",
      "gdT",
      "MAIT",
      "Th1-Like",
      "Th2-Activated",
      "Th2-Resting",
      "Th17-Resting",
      "Th17-Activated",
      "Th22",
      "Tfh-1",
      "Tfh-2",
      "Tph"
)
other_programs <- c(
      "Mito",
      "Poor-Quality",
      "IEG2",
      "IEG3",
      "Doublet-RBC",
      "Doublet-Platelet",
      "Doublet-Myeloid",
      "Doublet-Plasmablast",
      "Doublet-Bcell",
      "Doublet-Fibroblast"
)
all_programs <- c(activity_programs, identity_programs, other_programs)
# Check
sum(duplicated(all_programs))
setdiff(colnames(scores_tib), all_programs)

# Select relevant columns
metadata_programs_tib <- metadata_subset %>%
      select(
            patient_id,
            timepoint,
            celltype,
            cohort,
            TP53_status,
            all_of(activity_programs)
      )

# Checks
metadata_programs_tib %>%
      select(patient_id, cohort, TP53_status) %>%
      unique %>%
      count(cohort)
metadata_programs_tib$celltype %>% unique %>% length
length(activity_programs)

# Summarize to get one row per group (26 patients, 1 timepoint except P01 (2 timepoints), 10 celltypes except P18 (9 celltypes), 24 programs; 26*10*24+1*9*24 rows)
program_averages_tib <- metadata_programs_tib %>%
      pivot_longer(
            cols = all_of(activity_programs),
            names_to = "program",
            values_to = "usage"
      ) %>%
      group_by(
            patient_id,
            timepoint,
            cohort,
            celltype,
            program,
            TP53_status
      ) %>%
      summarize(mean_usage = mean(usage), .groups = "drop") %>%
      arrange(celltype) %>%
      mutate(celltype_patient = paste0(celltype, " (", patient_id, ")")) %>%
      mutate(
            celltype_patient = factor(
                  celltype_patient,
                  levels = unique(celltype_patient)
            )
      )

# Plot program usage per patient. The group_sizes is just to add lines in the plot below
group_sizes <- program_averages_tib %>%
      group_by(celltype, patient_id, cohort, TP53_status) %>%
      summarize(.groups = "drop") %>%
      group_by(celltype, cohort) %>%
      summarize(size = n(), .groups = "drop")
group_sizes <- group_sizes %>%
      group_by(cohort) %>%
      mutate(
            x_start = lag(cumsum(size), default = 0) + 1, # Start position for each group
            x_end = cumsum(size), # End position for each group
            x_center = (x_start + x_end) / 2 # Midpoint for labeling
      )

program_averages_tib %>%
      ggplot(aes(x = celltype_patient, y = program, fill = mean_usage)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      geom_segment(
            data = group_sizes,
            aes(
                  x = x_end + 0.5,
                  xend = x_end + 0.5,
                  y = 0.5,
                  yend = length(activity_programs) + 0.5
            ),
            color = "grey",
            linewidth = 0.5,
            inherit.aes = F
      ) +
      geom_hline(
            yintercept = 1:length(activity_programs) + 0.5,
            color = "grey"
      ) +
      facet_wrap(~cohort, scales = "free_x") +
      geom_text(
            data = group_sizes,
            aes(
                  x = x_center,
                  y = length(activity_programs) * 1.03,
                  label = celltype
            ),
            size = 3,
            angle = 45,
            hjust = 0,
            inherit.aes = F
      ) +
      scale_y_discrete(expand = expansion(mult = c(0, 0.22))) +
      theme_bw() +
      theme(
            panel.grid = element_blank(),
            aspect.ratio = 2,
            axis.text.x = element_text(angle = 45, hjust = 1, size = 5)
      )

ggsave("5.2.1_Program_usage_per_patient.png", width = 12, height = 8)

# The significance is hard to see in the plot above. Here's an alternative that shows it more clearly.
pdf("5.2.2_Program_usage_stats.pdf", width = 10, height = 8)

for (ct in unique(program_averages_tib$celltype)) {
      #ct <- "CD8 Effector"
      ct_averages_tib <- program_averages_tib %>% filter(celltype == ct)

      p_values <- ct_averages_tib %>%
            group_by(program) %>%
            summarize(
                  p_value = wilcox.test(
                        mean_usage ~ cohort,
                        data = cur_data()
                  )$p.value
            ) %>%
            ungroup() %>%
            mutate(p.adj = p.adjust(p_value, method = "BH")) %>%
            left_join(
                  group_by(ct_averages_tib, program) %>%
                        summarize(y_max = max(mean_usage) * 0.9)
            )

      p1 <- ct_averages_tib %>% #left_join(p_values)
            ggplot(aes(
                  x = cohort,
                  y = mean_usage,
                  color = cohort
            )) +
            geom_jitter(width = 0.2) +
            stat_compare_means(
                  method = "wilcox.test",
                  size = 3,
                  show.legend = F
            ) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
            scale_color_manual(values = c("#33CC0080", "#CE3D3280")) +
            facet_wrap(~program, scales = "free_y") +
            theme_bw() +
            ggtitle(ct) +
            theme(
                  panel.grid = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            geom_text(
                  data = p_values,
                  aes(
                        x = 1.5,
                        y = y_max,
                        label = paste0("p.adj = ", signif(p.adj, 3))
                  ),
                  inherit.aes = F,
                  size = 3
            )

      print(p1)
}

dev.off()

# Now show median program usage across patients, separated by cohort
program_averages_tib %>%
      group_by(cohort, celltype, program) %>%
      summarize(median_p = median(mean_usage), .groups = "drop") %>%
      ggplot(aes(x = celltype, y = program, fill = median_p)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      theme_minimal() +
      theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 2
      ) +
      facet_wrap(~cohort)

ggsave("5.2.3_Program_usage_medians.png", width = 12, height = 8)

# Note: in prior versions of this script, I also tried:
# 1. To subtract the medians and look at the difference between cohorts in one plot. This was misleading as small and insignificant changes were visually striking.
# 2. Making a heatmap to evaluate if samples could be clustered by T cell program usages. This clustering did not separate cohorts or TP53 mutation status
# I deleted these scripts on 250902 (should still be in Github history)
