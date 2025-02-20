# Ksenia Safina, 230908
# Edited by Peter van Galen, 250216
# To distinguish the host and donor cells, we have used the souporcell package in this analysis.
# Please visit https://github.com/wheaton5/souporcell to see how souporcell works.
# In this code, we have tested the quality output files of the souporcell.

# Load libraries
library(tidyverse)
library(fossil)
library(ggrepel)
library(ggpubr)
library(scales)

# Clear environment variables
rm(list=ls())

# Set working directory (local). For Nurefsan:
setwd("/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Souporcell/")
# For Peter:
#setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Souporcell/")

# Current patient to analyze
pt <- "P01"

# Create a file that combines reference and alternative matrices output of the souporcell in bash on Broad cluster:
#PT="P01"
#paste ref.mtx alt.mtx | sed 's/ /\t/g' > ${PT}.combined.tsv

# These files are saved in "5_Souporcell/AuxiliaryFiles" and not synced to GitHub due to their large size (see .gitignore). They are available upon request.
variants_df <- read.table(paste0("AuxiliaryFiles/", pt, ".combined.tsv"), skip = 3, sep = "\t")

# Name columns, remove variants without coverage, add genotype specifying if a variant is WT or MT
variants_df <- variants_df %>%
  dplyr::select(var = V1, cell = V2, ref = V3, alt = V6) %>%
  mutate(cov = ref + alt, genotype = ifelse(alt > 0, 0, 1)) %>%
  filter(cov > 0)

# Select barcodes and assignment from Soupourcell outputs
cells <- read.table(paste0("outputs/", pt, "_clusters.tsv"), header = T)
cells <- cells %>%
  dplyr::select(barcode, assignment)
cells$cell = 1:nrow(cells)

# Select cells with unambiguous souporcell assignment
cells$assignment %>% table()
cells = cells %>%
  filter(assignment %in% c(0,1))

# Add cell barcodes and souporcell assignments to variants dataframe. This is assuming combined.tsv and clusters.tsv have cells in the same order. Remove cells without unambiguous souporcell assignment.
variants_df <- variants_df %>% 
  left_join(cells) %>%
  drop_na()

# Extract variants with more than 100 cells, then filter for those variants
nice_variants <- variants_df %>% 
  group_by(var) %>%
  summarize(n = n()) %>%
  filter(n>=100) %>% pull(var)
nice_variants_df <- variants_df %>%
  filter(var %in% nice_variants)

# Show distribution of number of cells per variant
bind_rows(variants_df %>% count(var) %>% rename(`number of cells` = n) %>% mutate(Variants = "All"),
  nice_variants_df %>% count(var) %>% rename(`number of cells` = n) %>% mutate(Variants = "Nice")) %>%
  ggplot(aes(x = `number of cells`, fill = Variants)) +
  geom_histogram(position = "identity") +
  scale_x_log10() +
  ggtitle(paste0("Cells per variant (total = ", comma(length(unique(variants_df$cell))), " cells)")) +
  scale_fill_manual(values = c("lightblue", "lightgreen")) +
  theme_pubr()

# Variants per cell
bind_rows(variants_df %>% count(cell) %>%rename(`number of variants` = n) %>% mutate(Variants = "All"),
  nice_variants_df %>% count(cell) %>% rename(`number of variants` = n) %>% mutate(Variants = "Nice")) %>%
  ggplot(aes(x = `number of variants`, fill = Variants)) +
  geom_histogram(position = "identity") +
  scale_x_log10() +
  ggtitle(paste0("Variants per cell (total = ", comma(length(unique(variants_df$var))),
          "; ", comma(length(unique(nice_variants_df$var))), ")")) +
  scale_fill_manual(values = c("lightblue", "lightgreen")) +
  theme_pubr()

# For each variant, calculate the rand index to computes a similarity measure between the genotype and assignment. This can take a while.
rand_indices <- data.frame()
j <- 0
for (i in nice_variants) {
  j <- j+1
  cat(sprintf("\rPosition %s (%d/%d)", i, j, length(nice_variants)))
  genotypes <- nice_variants_df %>% filter(var == i) %>% pull(genotype)
  assignments <- nice_variants_df %>% filter(var == i) %>% pull(assignment) %>% as.numeric()
  r <- rand.index(genotypes, assignments)
  rand_indices <- rbind(rand_indices, data.frame(var = i , rand = r, ncells = length(genotypes)))
}

# Distribution of rand indices per variant
rand_indices %>%
  ggplot(aes(x = rand, y = ncells)) +
  geom_point() +
  theme_pubr()

# Add "geno_frac" to show the fraction of wt cells and mutated cells. Informative variants shouldn’t be close to 100%.
geno_frac <- nice_variants_df %>%
  group_by(var, genotype) %>%
  summarize(n=n()) %>%
  group_by(var) %>%
  summarize(max_geno_frac = max(n/sum(n)))

# Plot down below selects and label the variants that have more than 1000 cells and rand index >0.8 and geno fraction <0.8
rand_indices %>%
  left_join(geno_frac) %>%
  ggplot(aes(x = rand, y = ncells, color = max_geno_frac, label = var)) + 
  geom_point() +
  geom_text_repel(aes(label = ifelse(ncells>1000 & rand>0.8 & max_geno_frac<0.8,
                  as.character(var), "")), hjust = 0, vjust = 0, color = "black")+
  scale_color_distiller(palette = "Spectral") +
  theme_pubr()

# Select the most informative variants for to see how strong Souporcell results are 
best_variants <- rand_indices %>% left_join(geno_frac) %>%
  filter(rand > 0.8, ncells > 1000, max_geno_frac < 0.8) %>% pull(var)

# Arrange the x axes
ordered_barcodes <- cells %>% arrange(assignment) %>% pull(barcode)
nice_variants_df$barcode <- factor(nice_variants_df$barcode, levels = ordered_barcodes)

# Heatmap with top informative variants, cells are ordered by assignment column and colored by genotype (based on alt/ref counts)
p1 <- nice_variants_df %>%
  filter(var %in% best_variants) %>%
  mutate(genotype = as.character(genotype)) %>%
  ggplot(aes(x = barcode,y = as.character(var), fill = genotype)) + 
  geom_tile() +
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  labs(x = "Cell", y = "Variant", title = paste0(pt, " souporcell variants (",
    comma(length(unique(filter(nice_variants_df, var %in% best_variants)$cell))),
    " cells)")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        panel.grid = element_blank())

p1

# Save as a pdf file 
pdf(paste0("5.2_", pt, "_Heatmap.pdf"), width = 8, height = 4)
p1
dev.off()

# Save variants information
variants_df %>% select(barcode, var, ref, alt, genotype, assignment) %>% 
  filter(var %in% best_variants) %>%
  write_tsv(file = paste0("variants/", pt, "_variants.txt"))

