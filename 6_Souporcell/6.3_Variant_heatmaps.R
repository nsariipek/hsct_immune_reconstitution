# Ksenia Safina, Peter van Galen, updated 250302
# To distinguish the host and donor cells, we have used the souporcell package in this analysis (https://github.com/wheaton5/souporcell). Here, we assess the quality of variants detected by souporcell.

# Load libraries
library(tidyverse)
library(Seurat)
library(janitor)
library(fossil)
library(ggrepel)
library(ggpubr)
library(scales)
library(ComplexHeatmap)
library(cowplot)
library(grid)

# Clear environment variables
rm(list = ls())

# Favorite function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Set working directory (local). For Nurefsan:
setwd("/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Souporcell/")
# For Peter:
#setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/6_Souporcell/")

# Load Seurat data
seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")

# Load Souporcell calls
souporcell_assignments <- read_csv("6.2_Souporcell_assignments.csv.gz")

# Load cell type colors from 2.3_PvG-Colors.R
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

# Current patient(s) to analyze
pts <- sprintf("P%02d", 1:33)
#pt <- "P20" # for paper figure
#pt <- "P33" # faster example

# Open a connection to a stats file
file_conn <- file("6.3_stats_output.txt", open = "wt")
writeLines(
  paste(
    "Patient",
    "Total variants",
    "High coverage variants (max 1000)",
    "Most concordant variants",
    "Total cells",
    "Unambiguous cells",
    "Cells with most concordant variants",
    sep = "\t"
  ),
  file_conn
)
flush(file_conn)

# Loop over patients to generate heatmaps and other QC visualizations
for (pt in pts) {
  message(paste("Starting", pt, "analysis"))

  # For each patient, we made a file that combines reference and alternative matrices output of the souporcell in bash on Broad cluster in bash: `paste ref.mtx alt.mtx | sed 's/ /\t/g' > ${PT}.combined.tsv`
  # These files are saved in "5_Souporcell/AuxiliaryFiles" and not synced to GitHub due to their large size (see .gitignore). They are available upon request.
  variants_df <- read.table(paste0("AuxiliaryFiles/", pt, ".combined.tsv"), skip = 3, sep = "\t")

  # Name columns, remove variants without coverage, add genotype specifying if a variant is WT or MT
  variants_df <- variants_df %>%
    select(var = V1, cell = V2, ref = V3, alt = V6) %>%
    mutate(cov = ref + alt, genotype = ifelse(alt > 0, 0, 1)) %>%
    filter(cov > 0)
  message(paste("Total variants:", length(unique(variants_df$var))))

  # Select barcodes and assignment from Soupourcell outputs
  cells <- read.table(paste0("clusters/", pt, "_clusters.tsv"), header = T)
  cells <- cells %>%
    select(barcode, assignment)
  cells$cell = 1:nrow(cells)

  # Select cells with unambiguous souporcell assignment
  message(paste0(
    "Keeping cells with unambiguous souporcell call: ",
    nrow(filter(cells, assignment %in% c(0, 1))),
    "/",
    nrow(cells)
  ))
  unambiguous_cells <- cells %>%
    filter(assignment %in% c(0, 1))

  # Now, change Souporcell assignment (0 or 1) to origin (recipient or donor), or vice versa, as determined in 6.2_Genotype_donor-host_assignments.R
  zero_origin <- filter(souporcell_assignments, patient_id == pt, assignment == 0) %>% pull(origin) %>% unique
  one_origin <- filter(souporcell_assignments, patient_id == pt, assignment == 1) %>% pull(origin) %>% unique
  if(! "unknown" %in% c(zero_origin, one_origin)) {
    cells <- cells %>% mutate(assignment = gsub("0", substr(zero_origin, 1, 1), gsub("1", substr(one_origin, 1, 1), assignment))) %>%
      mutate(assignment = factor(assignment, levels = c("r", "r/d", "d", "d/r")))
    unambiguous_cells <- unambiguous_cells %>%
      mutate(assignment = gsub("0", zero_origin, gsub("1", one_origin, assignment))) %>%
      mutate(assignment = factor(assignment, levels = c("recipient", "donor")))
  }

  # Barplot of souporcell assignments
  p1 <- cells$assignment %>%
    tabyl() %>%
    rename("Number of cells" = n) %>%
    ggplot(aes(x = ., y = `Number of cells`)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = `Number of cells`), vjust = -0.3) +
    ylim(0, max(cells$assignment %>% tabyl() %>% pull(n)) * 1.1) +
    labs(
      x = "Souporcell assignment",
      title = paste0(
        "Cell assignments (",
        round(nrow(unambiguous_cells) / nrow(cells) * 100),
        "% unambiguous)"
      )
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank()
    )

  # Add cell barcodes and souporcell assignments to variants dataframe. This is assuming combined.tsv and clusters.tsv have cells in the same order. Remove cells without unambiguous souporcell assignment.
  variants_df <- variants_df %>%
    left_join(unambiguous_cells) %>%
    drop_na()

  # Select variants with coverage in at least 100 cells (no more than 1000 variants)
  high_coverage_variants <- variants_df %>%
    group_by(var) %>%
    summarize(n = n()) %>%
    filter(n > 100) %>%
    arrange(desc(n)) %>%
    head(1000) %>%
    pull(var)
  high_coverage_variants_df <- variants_df %>%
    filter(var %in% high_coverage_variants)

  # Show distribution of number of cells per variant
  p2 <- bind_rows(
    variants_df %>%
      count(var) %>%
      rename(`Number of cells` = n) %>%
      mutate(Variants = "All"),
    high_coverage_variants_df %>%
      count(var) %>%
      rename(`Number of cells` = n) %>%
      mutate(Variants = "High coverage")) %>%
    ggplot(aes(x = `Number of cells`, fill = Variants)) +
    geom_histogram(position = "identity") +
    scale_x_log10() +
    labs(y = "Number of variants") +
    ggtitle(paste0(
      "Cells per variant (total = ",
      comma(length(unique(variants_df$cell))),
      " cells)"
    )) +
    scale_fill_manual(values = c("lightblue", "darkgreen")) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank())

  # For each variant, calculate the rand index to computes a similarity measure between the genotype and assignment. This can take a while.
  rand_indices <- data.frame()
  j <- 0
  for (i in high_coverage_variants) {
    j <- j + 1
    cat(sprintf(
      "\rCalculating rand index for position %s (%d/%d)",
      i,
      j,
      length(high_coverage_variants)
    ))
    genotypes <- high_coverage_variants_df %>%
      filter(var == i) %>%
      pull(genotype)
    assignments <- high_coverage_variants_df %>%
      filter(var == i) %>%
      mutate(assignment = case_when(
        assignment %in% c("0", "donor") ~ 0,
        assignment %in% c("1", "recipient") ~ 1)) %>%
      pull(assignment)
    r <- rand.index(genotypes, assignments)
    rand_indices <- rbind(
      rand_indices,
      data.frame(var = i, rand = r, ncells = length(genotypes))
    )
  }

  # Select the most informative variants to see how strong Souporcell results are
  high_concordance_variants <- rand_indices %>% filter(rand > 0.9) %>% pull(var)
  high_concordance_variants_df <- filter(
    high_coverage_variants_df,
    var %in% high_concordance_variants
  )
  message(paste(
    "\nIdentified",
    length(high_concordance_variants),
    "concordant variants."
  ))

  # Plot rand index vs. number of cells
  p3 <- rand_indices %>%
    mutate(
      `Selected for heatmap` = factor(
        var %in% high_concordance_variants,
        levels = c("TRUE", "FALSE")
      )
    ) %>%
    ggplot(aes(x = rand, y = ncells, color = `Selected for heatmap`)) +
    geom_point() +
    scale_color_manual(values = c("TRUE" = "#00D68F", "FALSE" = "#D595A7FF")) +
    theme_bw() +
    labs(
      title = paste0(
        pt,
        " variants (",
        length(high_concordance_variants),
        "/",
        comma(nrow(rand_indices)),
        " are most concordant)"
      ),
      x = "Rand index",
      y = "Number of cells with coverage"
    ) +
    theme(panel.grid = element_blank(), aspect.ratio = 1)

  # Message how many cells will be missing from the heatmaps
  message(paste0(
    "Cells without variants will not be shown: ",
    sum(
      !unambiguous_cells$barcode %in%
        as.character(high_concordance_variants_df$barcode)
    ),
    "/",
    nrow(unambiguous_cells),
    " (",
    round(
      mean(
        !unambiguous_cells$barcode %in%
          as.character(high_concordance_variants_df$barcode)
      ) *
        100,
      2
    ),
    "%)"
  ))

  # Check if there are fewer than 100 cells with top variants (QC failed)
  if (
    sum(
      unambiguous_cells$barcode %in%
        as.character(high_concordance_variants_df$barcode)
    ) <
      100
  ) {
    message("Fewer than 100 cells...moving on to the next patient")

    # Save plots
    pdf(
      paste0("variant_heatmaps/", pt, "_QC_Failure.pdf"),
      width = 16,
      height = 3.5
    )
    print(plot_grid(p1, p2, p3, nrow = 1))
    dev.off()

    # Save stats
    writeLines(
      paste(
        pt,
        length(unique(variants_df$var)),
        length(high_coverage_variants),
        length(high_concordance_variants),
        nrow(cells),
        nrow(unambiguous_cells),
        sum(
          unambiguous_cells$barcode %in%
            as.character(high_concordance_variants_df$barcode)
        ),
        sep = "\t"
      ),
      file_conn
    )
    flush(file_conn)
    next
  }

  # Extract relevant metatadata from Seurat
  wrangled_metadata <- as_tibble(
    subset(seu, patient_id == pt)@meta.data,
    rownames = "barcode"
  ) %>%
    mutate(barcode = cutf(barcode, d = "_", f = 3)) %>%
    select(barcode, sample_status, celltype)

  # Join with variant level data
  heatmap_info_df <- high_concordance_variants_df %>%
    select(barcode, var, genotype, assignment) %>%
    left_join(wrangled_metadata)

  # To create complex heatmap with annotation bars, start with restructuring the data
  heatmap_wide_df <- heatmap_info_df %>%
    pivot_wider(
      id_cols = c(barcode, sample_status, assignment, celltype),
      names_from = var,
      values_from = genotype
    ) %>%
    as.data.frame() %>%
    mutate(
      sample_status = factor(
        sample_status,
        levels = c("pre-transplant", "remission", "relapse")
      )
    ) %>%
    sample_frac(1) %>%
    arrange(sample_status, assignment, celltype)

  # Heatmap annotation
  celltype_colors_current <- celltype_colors[factor(
    na.omit(unique(heatmap_wide_df$celltype)),
    levels = names(celltype_colors)
  )]
  col_anno <- columnAnnotation(
    sample_status = heatmap_wide_df$sample_status,
    souporcell_assignment = heatmap_wide_df$assignment,
    celltype = heatmap_wide_df$celltype,
    col = list(
      sample_status = c(
        "pre-transplant" = "#A3BFD9",
        "remission" = "#F6E06E",
        "relapse" = "#8B0000"
      ),
      souporcell_assignment = c("0" = "#3B1B53", "1" = "#F0E685",
        "donor" = "#4B3140", recipient = "#E4C9B0"),
      celltype = celltype_colors_current
    ),
    annotation_legend_param = list(
      celltype = list(
        ncol = 2
      )
    )
  )

  # Convert numeric matrix to character matrix with "0" and "1"
  heatmap_mat <- heatmap_wide_df %>%
    column_to_rownames("barcode") %>%
    select(-c("sample_status", "assignment", "celltype")) %>%
    t()

  # Generate heatmap
  h2 <- Heatmap(
    heatmap_mat,
    col = c("0" = "#56B4E9", "1" = "#E69F00"),
    na_col = "white",
    cluster_rows = F,
    cluster_columns = F,
    use_raster = ncol(heatmap_mat) * nrow(heatmap_mat) > 100000,
    raster_quality = 5,
    bottom_annotation = col_anno,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 0),
    column_title = paste0(
      comma(ncol(heatmap_mat)),
      "/",
      comma(length(unique(unambiguous_cells$barcode))),
      " cells from ",
      pt,
      " (",
      round(
        ncol(heatmap_mat) / length(unique(unambiguous_cells$barcode)) * 100
      ),
      "%)"
    ),
    row_title = paste(nrow(heatmap_mat), "variants"),
    heatmap_legend_param = list(title = "genotype"),
    border = T
  )

  # Create the top row plot first
  top_row <- plot_grid(p1, p2, p3, ncol = 3)
  top_grob <- ggplotGrob(top_row)

  # Save the combined plot
  pdf(paste0("variant_heatmaps/", pt, "_Heatmap.pdf"), width = 16, height = 10)
  # Plot the heatmap, leaving space for the top row
  draw(h2, padding = unit(c(0.2, 0.2, 4, 0.2), "inch"))

  # Add the top row above the heatmap
  pushViewport(viewport(x = -0.02, y = 0.98, height = 0.35, just = c(0, 1)))
  # Add the top row as a grob of a ggplot object
  grid.draw(top_grob)

  dev.off()

  # Write statistics to the output files
  writeLines(
    paste(
      pt,
      length(unique(variants_df$var)),
      length(high_coverage_variants),
      length(high_concordance_variants),
      nrow(cells),
      nrow(unambiguous_cells),
      sum(
        unambiguous_cells$barcode %in%
          as.character(high_concordance_variants_df$barcode)
      ),
      sep = "\t"
    ),
    file_conn
  )
  flush(file_conn)

  message(paste("Done with", pt))

  # Trigger garbage collection
  gc()

}

# Close the file connection
close(file_conn)

# Looking at the heatmaps, there's a population of cells in Patient 24 that was likely misassigned by Souporcell. Patient 20 and 30 look particularly good.