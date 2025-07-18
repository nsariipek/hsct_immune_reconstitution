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
#VM: setwd("~/TP53_ImmuneEscape/3_DGE/")
#Local: setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/3_DGE")

# Delete environment variables & load favorite function
rm(list = ls())
cutf <- function(x, f = 1, d = "/") {
      sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")
t_celltypes <- c(
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
seu_T <- subset(seu, celltype %in% t_celltypes)

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

# Add T cell UMAP coordinates generated in 4_Trajectories/4.2_TCR_Diversity_UMAP.R & check visually
t_coordinates <- read.table("../4_Trajectories/4.2_UMAP-embeddings.csv", header = T, sep = ",", row.names = 1)
seu_T[["umapT"]] <- CreateDimReducObject(embeddings = as.matrix(t_coordinates), key = "umapT_", assay = "RNA")
DimPlot(seu_T, reduction = "umapT", shuffle = T, group.by = "TCAT_Multinomial_Label", cols = celltype_colors) +
      theme_bw() +
      theme(aspect.ratio = 1, panel.grid = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 3)))


### PETER LEFT OFF HERE


# Combine TCAT scores and program usage results with Seurat metadata
usage_tib <- read_tsv("3.1_starCAT/starCAT.scores.txt.gz") %>%
      rename("cell" = "...1")
scores_tib <- read_tsv("3.1_starCAT/starCAT.rf_usage_normalized.txt.gz") %>%
      rename("cell" = "...1")
metadata_tib <- as_tibble(seu_T@meta.data, rownames = "cell")
metadata_tib <- left_join(metadata_tib, scores_tib)
metadata_tib <- left_join(metadata_tib, usage_tib)

# Filter data for relevant time points and sample status
metadata_100d_tib <- metadata_tib %>%
      filter(
            timepoint %in% c(3, 5, 6),
            sample_status == "remission"
      )

# Quick look at antigen-specific activation (ASA) scores...Is the difference between cohorts driven by TP53?
p1 <- metadata_100d_tib %>%
      filter(celltype == "CD8 Effector Memory 2") %>%
      ggplot(aes(
            x = cohort,
            y = `ASA`,
            color = TP53_status == "MUT"
      )) +
      geom_violin(scale = "width", fill = NA, draw_quantiles = 0.5) +
      stat_compare_means(method = "wilcox.test", show.legend = F) + # This may not be the best statistical test
      theme_bw() +
      #nurefsan's edition
      theme(panel.grid = element_blank())
# Save the plot cohort# Save the plot
pdf("3.2.1_ASA_Score.pdf", width = 8, height = 6)
p1
dev.off()

# Subset for analysis of activity programs from the paper (Figure 2C, Table S1)
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
metadata_programs_tib <- metadata_100d_tib %>%
      select(
            patient_id,
            timepoint,
            celltype,
            cohort,
            TP53_status,
            all_of(activity_programs)
      )

# Summarize to get mean program usage per patient per cell type
metadata_programs_tib %>%
      select(patient_id, cohort, TP53_status) %>%
      unique %>%
      count(cohort)
metadata_programs_tib$celltype %>% unique %>% length
length(activity_programs)
# Here we expect a tibble with nrow is the number of patients * cell types & programs
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
            x_center = (x_start + x_end) / 2
      ) # Midpoint for labeling

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

ggsave("3.2.2_PerPatient.png", width = 12, height = 8)

# The significance is hard to see in the plot above. Here's an alternative that shows it more clearly.
pdf("3.2.3_Test_plots.pdf", width = 10, height = 8)

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

# Now show the median across patients. Of note, this is misleading, as visually striking differences may in fact be driven by 1-2 patients
program_averages_tib %>%
      group_by(cohort, celltype, program) %>%
      summarize(median_p = median(mean_usage), .groups = "drop") %>% # Mean program usage per cohort (3)
      ggplot(aes(x = celltype, y = program, fill = median_p)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      theme_minimal() +
      theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 2
      ) +
      facet_wrap(~cohort)

# Subtract to identify look at differences
non_relapsed_tib <- program_averages_tib %>%
      filter(cohort == "long-term-remission") %>%
      group_by(cohort, celltype, program) %>%
      summarize(median_p = median(mean_usage), .groups = "drop")

relapsed_tib <- program_averages_tib %>%
      filter(cohort == "relapse") %>%
      group_by(cohort, celltype, program) %>%
      summarize(median_p = median(mean_usage), .groups = "drop")

join_tib <- full_join(
      non_relapsed_tib,
      relapsed_tib,
      by = c("celltype", "program"),
      suffix = c(".non-relapsed", ".relapsed")
)
join_tib <- join_tib %>%
      mutate(diff = median_p.relapsed - `median_p.non-relapsed`)

join_tib %>%
      ggplot(aes(x = celltype, y = program, fill = diff)) +
      geom_tile() +
      scale_fill_gradientn(
            colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 2)
ggsave("3.2.4_Difference.png", width = 12, height = 8)

# Now make a heatmap like Figure 4F in Van Galen ... Bernstein 2019

# Convert data to matrix
# 9 cell types, 24 programs --> 216 rows
df <- program_averages_tib %>%
      mutate(celltype_program = paste0(celltype, ", ", program)) %>%
      mutate(
            cohort_patient = paste(
                  cohort,
                  TP53_status,
                  patient_id,
                  ".",
                  timepoint,
                  "M"
            )
      ) %>%
      select(celltype_program, cohort_patient, mean_usage) %>%
      pivot_wider(names_from = cohort_patient, values_from = mean_usage) %>%
      as.data.frame() %>%
      column_to_rownames("celltype_program")
df[1:3, 1:3]

# Side annotation
T_cell_types <- intersect(levels(seu_T$celltype), seu_T$celltype)
T_cell_types <- factor(T_cell_types, levels = T_cell_types)
row_anno <- rowAnnotation(
      celltype = factor(
            cutf(rownames(df), d = ", ", f = 1),
            levels = levels(T_cell_types)
      ),
      program = cutf(rownames(df), d = ", ", f = 2),
      col = list(celltype = celltype_colors)
)

# Bottom annotation
tokens <- strsplit(colnames(df), " ")
cohort_vec <- sapply(tokens, `[`, 1)
TP53_vec <- sapply(tokens, `[`, 2)

# Reorder df by cohort
cohort_factor <- factor(
      cohort_vec,
      levels = c("long-term-remission", "relapse")
)
column_order <- order(cohort_factor)
df_ordered <- df[, column_order]
# Apply same order to annotation variables
cohort_vec_ordered <- cohort_vec[column_order]
TP53_vec_ordered <- TP53_vec[column_order]

cohort_colors <- c("long-term-remission" = "#546fb5FF", "relapse" = "#e54c35FF")
TP53_colors <- c("WT" = "#6BD76BFF", "MUT" = "#FF1463FF")
col_anno <- HeatmapAnnotation(
      Cohort = factor(
            cohort_vec_ordered,
            levels = c("long-term-remission", "relapse")
      ),
      TP53_status = factor(TP53_vec_ordered, levels = c("WT", "MUT")),
      col = list(
            Cohort = cohort_colors,
            TP53_status = TP53_colors
      )
)

# Generate heatmap with clustering
h1 <- Heatmap(
      as.matrix(df_ordered),
      col = colorRamp2(
            c(
                  min(as.matrix(df_ordered), na.rm = TRUE),
                  max(as.matrix(df_ordered), na.rm = TRUE)
            ),
            c("white", "red")
      ),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      left_annotation = row_anno,
      bottom_annotation = col_anno,
      row_names_gp = gpar(fontsize = 6)
)

h1

png("3.2.5_ClusteredHeatmap.png", height = 1500, width = 800)
draw(h1)
dev.off()
