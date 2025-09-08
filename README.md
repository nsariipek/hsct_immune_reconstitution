# Immune Reconstitution After Hematopoietic Stem Cell Transplantation

This repository contains single-cell analysis from 33 AML/MDS patients over the course of hematopoietic stem cell transplantation (HSCT). We performed paired 10x Chromium Single Cell 5' RNA-sequencing and V(D)J Enrichment for 48 unique samples and 65 libraries in total.

## 00_PreProcessing

To generate count matrices, we first created configuration files according to [10x instructions](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-5p-multi) (1_CreateConfigFiles.sh). We then submitted Cell Ranger jobs with standard GRCh38 and VDJ references (2_RunCellRangerMulti.sh, 3_Cellranger.sh).

## 01_Seurat

We created Seurat object from gene expression data from all samples and added sample information such as patient identification and time point. We assessed cell numbers, visualized quality control metrics, and filtered for high-quality cells (1.1_CreateSeuratObject.R).

## 02_Annotate-predict

For cell annotation, we ran cell type prediction using [BoneMarrowMap](https://github.com/andygxzeng/BoneMarrowMap). Each of the samples was mapped onto the BoneMarrowMap reference, annotations were merged, and some broad cell type labels were replaced by more granular cell type labels (2.1_Project_BoneMarrowMap_annotation.R). In 2.2_Complete_Seurat_object.R, we integrated all different modalities (including annotations, TCR calls, Souporcell and Numbat results) into one final Seurat object saved as `250528_Seurat_complete.rds`. This file is the basis of most downstream analyses and will be available on Figshare upon publication. In 2.3_MarkerGeneHeatmap.R, we created a heatmap of cell type markers.

## 03_Cell_Proportions

This folder contains scripts for cell type proportion analysis, including stacked barplots, T cell proportions, and the CD4/CD8 ratio across patients.

## 04_Trajectories

We performed trajectory analysis of T cell subsets using [Monocle3](https://github.com/cole-trapnell-lab/monocle3). The complete Seurat object was subsetted to non-proliferating CD4 or CD8 T cells, followed by Seurat dimensionality reduction (PCA and Harmony), clustering, and UMAP visualization. The resulting Seurat object is converted to a Monocle3 cell_data_set for trajectory inference, with naive T cell clusters set as the root. Pseudotime values are then calculated, visualized (UMAPs, boxplots, density plots), and compared between remission and relapse cohorts. For pseudotime analysis of malignant cells, see 09_Numbat.

## 05_DGE

Differential gene expression (DGE) analysis was performed for T cells using [starCAT](https://github.com/immunogenomics/starCAT) (5.1_starCAT.sh, 5.1_TCAT.R, and 5.2_TCAT_Programs.R). DGE on recipient hematopoietic stem and progenitor cells (HSPCs) vs. their donor counterparts was done on pseudobulked cells using DESeq2 and apeglm packages. This was followed by fast gene set enrichment (FGSEA) analysis on genes sorted by log2 fold change (5.3_DGE_HSPCs.R and 5.4_GSEA_HSPCs.R). DGE on relapse vs. pre-transplant samples was similarly done on pseudobulked cells using DESeq2 and apeglm packages. Thereafter, we again performed FGSEA analysis on genes sorted by log2 fold change (5.5_DGE_over-time.R and 5.6_GSEA_over-time.R). We separately checked expression of MHC-II related genes (5.7_MHCII_over-time.R).

## 06_TCR_Diversity

We combined gene expression data (Seurat object metadata) with TCR data by loading the Cell Ranger Multi output /filtered_contig_annotations.csv files for each sample (6.1_Integrate_TCR_Seurat.R). We used the [scRepertoire](https://github.com/BorchLab/scRepertoire) for subsequent analysis and adding the `CTStrict` column to the Seurat object for downstream TCR analyses. In 6.4_TCR_Diversity_3M.R, we compared TCR diversity between long-term remission and relapse cohorts, focusing specifically on remission samples \~3 months after transplant, downsampling T cells and calculating their diversity using the Inverse Simpson Index (DiversityFunctions.R). Other scripts in this folder visualize clonotype dominance (6.2_TCR_Diversity_UMAPs, 6.3_Clonotype_dominance_3M.R). TCR diversity over time is visualized using 6.5_Longitudinal_diversity.R

## 07_Neoantigen

We first subset the complete Seurat object to T cells. We obtained antigen recognition signatures from published datasets (Rosenberg Lab, [Lowery 2022](https://www.ncbi.nlm.nih.gov/pubmed/35113651), separate signatures for CD4+ and CD8+ T cells; [Yossef 2023](https://www.ncbi.nlm.nih.gov/pubmed/38039963) for CD8+ T cells), scored single cells, and calculated the median score per T cell clonotype and per patient. We used the neoantigen-specific activation (ASA) scores from [TCAT](https://www.ncbi.nlm.nih.gov/pubmed/40903640) for CD4+ and CD8+ T cells in a similar manner. Data were visualized as sina/violin plots ordered from high to low median scores with tiles indicating cohort and TP53 mutation status, as well as box plots comparing the medians between patients with TP53-wildtype vs. TP53-mutated AML/MDS.

## 08_Souporcell

We distinguished recipient and donor cells using the [Souporcell](https://github.com/wheaton5/souporcell) package. Barcodes were first extracted from the Seurat object for each patient; in longitudinal cases, barcodes from multiple timepoints were merged and duplicates removed (8.1_Extract_barcodes_for_Souporcell.R). Souporcell was then run on a Google virtual machine and output (clusters) extracted. We applied rule based selection to assign genotypes 0 and 1 to recipient and donor, and added the assignments to the Seurat object (8.2_Add_Souporcell_assignments_to_metadata.R). Top variants were selected and visualized in 8.3_Variant_heatmaps.R. The proportions of recipient and donor cells in different patients, time points, and cell types were visualized in 8.4_Souporcell_plots.R, including HSPC chimerism plots.

## 09_Numbat

We further explored our dataset using [Numbat](https://github.com/kharchenkolab/numbat) to identify malignant (tumor) cells, focusing only on recipient cells as inferred by Souporcell. Recipient barcodes and corresponding counts were extracted (9.1_Extracting_barcodes.R, 9.2_Extracting_countdata.R) and used to run the analysis and a computing cluster. The resulting CNV calls and plotting data were transferred back and visualized in R (9.3_Visualize_Numbat_results.R), followed by integration into the Seurat object using the output file PX_clone_post_2.tsv (9.4_Integrate_Numbat_Seurat.R). Confidence levels of the CNV calls were evaluated using metrics provided by Numbat (9.5_confidence_level.R). The malignant cell proportion across cell types was visualized in 9.6_Malignant_proportion_heatmap.R. In 9.7_Other_plots.R, we generate UMAPs and determine that remission HSPCs harbor CNVs. Finally, BoneMarrowMap-derived pseudotime of malignant cells was compared between pre-transplant and relapse samples in 9.8_Malignant_cell_pseudotime.R.

## 10_Miscellaneous

Additional code to produce specific figure panels: pie chart of cohort sizes, swimmer plot, and clinical chimerism box plot.