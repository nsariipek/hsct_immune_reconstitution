# Peter van Galen, 250206
# Assess TCAT programs in T cells

library(Seurat)
library(tidyverse)
library(scattermore)
library(ggpubr)
library(ggforce)
library(ComplexHeatmap)
library(circlize)

# Set working directory
setwd("~/TP53_ImmuneEscape/3_DGE/")

# Delete environment variables & load favorite function
rm(list=ls())
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load data
seu <- readRDS("~/250128_seurat_annotated_final.rds")

# Load colors from 2.3_PvG-Colors.R
celltype_colors_df <- read.table("../celltype_colors.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)

# Subset Seurat object to T cells
T_cell_types <- c("CD4 Naïve", "CD4 Effector Memory", "CD4 Memory", "Treg", "CD8 Naïve", "CD8 Effector", "CD8 Memory", "CD8 Exhausted", "γδ T")
T_cell_types <- factor(T_cell_types, levels = T_cell_types)
seu_T <- subset(seu, celltype %in% T_cell_types)

# Check visually
seu_T@meta.data %>% sample_frac(1) %>%
          ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = celltype)) +
          geom_scattermore(pointsize = 8, pixels = c(4096, 4096)) +
          scale_color_manual(values = celltype_colors) +
          theme_bw() +
          theme(aspect.ratio = 1,
                panel.grid = element_blank()) +
          guides(color = guide_legend(override.aes = list(size = 3)))

# Combine TCAT scores and program usage results with Seurat metadata
usage_tib <- read_tsv("AuxiliaryFiles/results.rf_usage_normalized.txt") %>% rename("cell" = "...1")
scores_tib <- read_tsv("AuxiliaryFiles/results.scores.txt") %>% rename("cell" = "...1")
metadata_tib <- as_tibble(seu_T@meta.data, rownames = "cell")
metadata_tib <- left_join(metadata_tib, scores_tib)
metadata_tib <- left_join(metadata_tib, usage_tib)

# Optional:
TP53_mut <- c("P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P12", "P14", "P17")
#metadata_tib <- filter(metadata_tib, patient_id %in% TP53_mut)
#metadata_tib <- filter(metadata_tib, ! patient_id %in% TP53_mut)

# Filter data for relevant time points etc.
metadata_100d_tib <- metadata_tib %>%
      filter(timepoint %in% c("3","5","6"),  # 100 days after transplant
             sample_status == "remission") %>%
      mutate(cohort_binary = substr(cohort, 3, nchar(cohort)))

# Quick look at antigen-specific activation (ASA) scores...Is the difference between cohorts driven by TP53?
metadata_100d_tib %>% filter(celltype == "CD8 Effector") %>%
      ggplot(aes(x = cohort_binary, y = `ASA`, color = patient_id %in% TP53_mut)) +
      geom_violin(scale = "width", fill = NA, draw_quantiles = 0.5) +
      theme_bw()

# Subset for analysis of activity programs from the paper (Figure 2C, Table S1)
programs <- c("BCL2/FAM13A", "TIMD4/TIM3", "ICOS/CD38", "CD172a/MERTK", "SOX4/TOX2", "Heatshock",
              "OX40/EBI3", "CD40LG/TXNIP", "Exhaustion", "CellCycle-G2M", "CellCycle-S",
              "Cytoskeleton", "Cytotoxic", "ISG", "CellCycle-Late-S", "HLA", "NME1/FABP5",
              "Translation", "IEG", "IL10/IL19", "RGCC/MYADM", "Metallothionein",
              "Multi-Cytokine", "CTLA4/CD38")
#TEST: metadata_100d_tib$celltype <- factor(metadata_100d_tib$Multinomial_Label, levels = c("CD4_Naive", "CD4_EM", "CD4_CM", "Treg", "CD8_Naive", "CD8_EM", "CD8_CM", "CD8_TEMRA", "gdT", "MAIT"))
metadata_programs_tib <- metadata_100d_tib %>%
      select(patient_id, celltype, cohort_binary, all_of(programs))

# Summarize to get mean program usage per patient per cell type
metadata_programs_tib$patient_id %>% unique %>% length # 28 patients
metadata_programs_tib %>% select(patient_id, cohort_binary) %>% unique %>% count(cohort_binary) # 18 non-relapsed, 10 relapsed
metadata_programs_tib$celltype %>% unique %>% length # 9 cell types
length(programs) # 24 programs
# So we expect 28*9*24=6,048 rows from this:
program_averages_tib <- metadata_programs_tib %>%
      pivot_longer(cols = programs, names_to = "program", values_to = "usage") %>%
      group_by(patient_id, cohort_binary, celltype, program) %>%
      summarize(mean_usage = mean(usage), .groups = "drop") %>%
      arrange(celltype) %>%
      mutate(celltype_patient = paste0(celltype, " (", patient_id, ")")) %>%
      mutate(celltype_patient = factor(celltype_patient, levels = unique(celltype_patient))) 

# Plot program usage per patient. The group_sizes is just to add lines in the plot below
group_sizes <- program_averages_tib %>% group_by(celltype, patient_id, cohort_binary) %>% summarize(.groups = "drop") %>%
      group_by(celltype, cohort_binary) %>% summarize(size = n(), .groups = "drop")
group_sizes <- group_sizes %>% group_by(cohort_binary) %>%
  mutate(x_start = lag(cumsum(size), default = 0) + 1,  # Start position for each group
         x_end = cumsum(size),  # End position for each group
         x_center = (x_start + x_end) / 2)  # Midpoint for labeling

program_averages_tib %>%
      ggplot(aes(x = celltype_patient, y = program, fill = mean_usage)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      geom_segment(data = group_sizes, aes(x = x_end + 0.5, xend = x_end + 0.5, y = 0.5, yend = 24.5),
                   color = "grey", linewidth = 0.5, inherit.aes = F) +
      geom_hline(yintercept = 1:24+0.5, color = "grey") +
      facet_wrap(~ cohort_binary, scales = "free_x") +
      geom_text(data = group_sizes, aes(x = x_center, y = 24.7, label = celltype), size = 3, angle = 45, hjust = 0, inherit.aes = F) +
      scale_y_discrete(expand = expansion(mult = c(0, 0.22))) +
      theme_bw() +
      theme(panel.grid = element_blank(), aspect.ratio = 2,
            axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
      
#ggsave("3.3.1_PerPatient.png", width = 12, height = 8)

# The significance is hard to see in the plot above. Here's an alternative that shows it more clearly.
pdf("3.3.2_Test_plots.pdf", width = 10, height = 8)

for ( ct in unique(program_averages_tib$celltype) ) {
      #ct <- "CD8 Effector"
      ct_averages_tib <- program_averages_tib %>% filter(celltype == ct)

      p_values <- ct_averages_tib %>% group_by(program) %>%
            summarize(p_value = wilcox.test(mean_usage ~ cohort_binary, data = cur_data())$p.value) %>% ungroup() %>%
            mutate(p.adj = p.adjust(p_value, method = "BH")) %>%
            left_join(group_by(ct_averages_tib, program) %>% summarize(y_max = max(mean_usage) * 0.9))
      
      p1 <- ct_averages_tib %>% #left_join(p_values)
            ggplot(aes(x = cohort_binary, y = mean_usage, color = cohort_binary)) +
            geom_jitter(width = 0.2) +
            stat_compare_means(method = "wilcox.test", size = 3, show.legend = F) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
            scale_color_manual(values = c("#33CC0080", "#CE3D3280")) +
            facet_wrap(~ program, scales = "free_y") +
            theme_bw() +
            ggtitle(ct) +
            theme(panel.grid = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1)) +
            geom_text(data = p_values, aes(x = 1.5, y = y_max,
                  label = paste0("p.adj = ", signif(p.adj, 3))), inherit.aes = F, size = 3) 
      
      print(p1)
}

dev.off()

# Now show the median across patients
program_averages_tib %>%
      group_by(cohort_binary, celltype, program) %>%
      summarize(median_p = median(mean_usage), .groups = "drop") %>% # Mean program usage per cohort (3)
            ggplot(aes(x = celltype, y = program, fill = median_p)) +
            geom_tile() +
            scale_fill_gradient(low = "white", high = "red") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 2) +
            facet_wrap(~ cohort_binary)

#ggsave("3.3.3_Program-usage-cohorts.png", width = 8, height = 8)

# Subtract to identify look at differences
non_relapsed_tib <- program_averages_tib %>% filter(cohort_binary == "Non-relapsed") %>%
      group_by(cohort_binary, celltype, program) %>%
      summarize(median_p = median(mean_usage), .groups = "drop")

relapsed_tib <- program_averages_tib %>% filter(cohort_binary == "Relapsed") %>%
      group_by(cohort_binary, celltype, program) %>%
      summarize(median_p = median(mean_usage), .groups = "drop")

join_tib <- full_join(non_relapsed_tib, relapsed_tib, by = c("celltype", "program"),
                      suffix = c(".non-relapsed", ".relapsed"))
join_tib <- join_tib %>% mutate(diff = median_p.relapsed - `median_p.non-relapsed`)

join_tib %>%
      ggplot(aes(x = celltype, y = program, fill = diff)) +
      geom_tile() +
      scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 2)

#ggsave("3.3.4_Program-usage-diff.png", width = 5, height = 8)




# UNDER CONSTRUCTION ↓↓↓

# Now make a heatmap like Figure 4F in Van Galen ... Bernstein 2019

# Convert data to matrix
# 9 cell types, 24 programs --> 216 rows

df <- program_averages_tib %>% mutate(celltype_program = paste0(celltype, ", ", program)) %>%
      mutate(cohort_patient = paste(cohort_binary, patient_id)) %>%
      select(celltype_program, cohort_patient, mean_usage) %>%
      pivot_wider(names_from = cohort_patient, values_from = mean_usage) %>%
      as.data.frame() %>%
      column_to_rownames("celltype_program")
df[1:3,1:3]

# Side annotation
row_anno <- rowAnnotation(celltype = factor(cutf(rownames(df), d = ", ", f = 1), levels = levels(T_cell_types)),
                          program = cutf(rownames(df), d = ", ", f = 2),
                          col = list(celltype = celltype_colors))
col_anno <- columnAnnotation(TP53_mut = grepl(paste(TP53_mut, collapse = "|"), colnames(df)),
                             col = list(TP53_mut = c(`TRUE` = "#FF1463FF", `FALSE` = "#6BD76BFF")))

# Generate heatmap with clustering
h1 <- Heatmap(as.matrix(df),
        col = colorRamp2(c(min(as.matrix(df), na.rm = TRUE), max(as.matrix(df), na.rm = TRUE)), c("white", "red")),
        cluster_rows = FALSE,
        left_annotation = row_anno,
        bottom_annotation = col_anno,
        row_names_gp = gpar(fontsize = 6))

h1

png("3.3.5_ClusteredHeatmap.png", height = 1500, width = 800)
draw(h1)
dev.off()