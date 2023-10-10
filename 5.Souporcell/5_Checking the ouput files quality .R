#Ksenia Safina, 230908
#To distinguish the host and donor cells, we have used the souporcell package in this analysis.
#Please visit https://github.com/wheaton5/souporcell to see how souporcell works.
#In this code, we have tested the quality output files of the souporcell.

#Load libraries
library(tidyverse)
library(fossil)
library(ggrepel)


# Start with a clean slate
rm(list=ls())

setwd("/Users/dz855/Downloads/")

#Create a file that combines reference and alternative matrices output of the souporcell.
# The way to create this file in bash is
# paste ref.mtx alt.mtx | sed 's/ /\t/g' > Pt09.combined.tsv

# t = read.table("ptbm.combined.mtx2", skip=3, sep = '\t')

t = read.table("souporcell_Pt9_combined.tsv", skip=3, sep = '\t')


#muatate the coverage name, add a genotype column specifying if a variant is WT or MT
t = t %>%
  dplyr::select(pos=V1, cell = V2, ref = V3, alt = V6) %>%
  mutate(cov = ref + alt, genotype = ifelse(alt > 0, 0,1)) %>%
  filter(cov > 0)

#add the output of the souporcell, and select the barcode and assignment.
cells = read.table("souporcell_Pt9_clusters.tsv", header = T)
cells = cells %>%
  dplyr::select(barcode, assignment)


cells$cell = 1:nrow(cells)

#select the cells only assigned correctly by souporcell
cells = cells %>%
  filter(assignment %in% c(0,1))

t = t %>% 
  left_join(cells) %>%
  drop_na()

#call variants according to the unique position of each variance and filter the ones that have less than 100 cells.

variants = unique(t$pos)

variants = t %>% 
  group_by(pos) %>%
  summarize(n = n()) %>%
  filter(n>=100) %>% pull(pos)

t_f = t %>%
  filter(pos %in% variants)

#For each variance, calculate the rand index, which computes a similarity measure between 2 clusterings, in this case, genotype and assignment.

results = data.frame()

for (i in variants) {
  print(i)
  genotypes = t_f %>% filter(pos == i) %>% pull(genotype)
  assignments = t_f %>% filter(pos == i) %>% pull(assignment) %>% as.numeric()
  r = rand.index(genotypes, assignments)
  results = rbind(results, data.frame(var= i ,rand = r, ncells = length(genotypes)))
}

t %>% 
  group_by(pos) %>%
  summarize(n = n()) %>%
  ggplot(aes(x=n)) + geom_histogram() + scale_x_log10()

results %>%
  ggplot(aes(x=rand,y=ncells)) + geom_point()

#add "geno_frac" showing the fraction of wt cells fraction of mutated cells. Informative variants shouldn’t be close to 100%.

geno_frac = t_f %>%
  group_by(pos, genotype) %>%
  summarize(n=n()) %>%
  group_by(pos) %>%
  summarize(max_geno_frac = max(n/sum(n)))

results %>% 
  left_join(geno_frac, by = c("var" = "pos")) %>%
  ggplot(aes(x=rand,y=ncells,color=max_geno_frac, label = var, size =4)) + 
  geom_point()+ 
  theme_pubr() +
  geom_label_repel(aes(label=ifelse(ncells>1000 & rand>0.8 & max_geno_frac<0.8, as.character(var),'')),hjust=0,vjust=0, label.size = 1)+
  scale_color_distiller(palette = "Spectral")



relevant_variants = results %>% left_join(geno_frac, by = c("var" = "pos")) %>% filter(rand > 0.8, ncells > 1000, max_geno_frac<0.8) %>% pull(var)

relevant_variants = results %>% left_join(geno_frac, by = c("var" = "pos")) %>% filter(rand > 0.8, ncells > 200, max_geno_frac<0.8) %>% pull(var)


ordered_barcodes = cells %>% arrange(assignment) %>% pull(barcode)

t_f$barcode = factor(t_f$barcode, levels = ordered_barcodes)


#Heatmap with top informative variants, cells are ordered by assignment column and colored by genotype (based on alt/ref counts)
t_f %>%
  filter(pos %in% relevant_variants) %>%
  mutate(genotype = as.character(genotype)) %>%
  dplyr::select(pos, genotype, barcode, assignment) %>% 
  ggplot(aes(x=barcode,y=as.character(pos),fill=genotype)) + 
  geom_tile() +
  theme(axis.text.x = element_blank())


df <- subset(t_f, select = -c(alt,ref))

rownames(df)  <- df$unique(barcode)      

#see the number of cells that are plotted(unique barcodes)
  t_f %>% 
  filter(pos %in% relevant_variants) %>%
  mutate(genotype = as.character(genotype)) %>%
  dplyr::select(pos, genotype, barcode, assignment) %>% pull(barcode) %>% unique() %>% length()


t %>% filter(pos %in% relevant_variants) %>% View()










