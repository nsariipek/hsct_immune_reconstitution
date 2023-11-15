# TP53_ImmuneEscape

This repository contains single-cell analysis on bone marrow samples from 12 TP53 mutated-AML patients longitudinally during the course of stem cell transplantation. We performed paired 10x Chromium Single Cell 5' and V(D)J Enrichment for 27 unique samples 44 libraries in total. Single-cell transcriptomes were aligned to hg38 using custom scripts Cell Ranger (10x data). 

Load gene expression data from samples, integrate, reduce dimensionality, cluster (1.Seurat)

Annotate cell clusters using canonical marker genes (2.Annotate)

Add single-cell genotyping metadata to the Seurat objects containing gene expression and cell type annotations (3.DGE)

                                                                                          (4.GSEA)

Distuingish cells as donor and host cells using Souporcell package (Heaton,2019) (5.Souporcell)

Analysis of TCR diversity from TCR-Seq(6.TCR Diversity)
