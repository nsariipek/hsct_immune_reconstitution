# Nurefsan Sariipek, 240823
# Analyzing samples from 1801 project from Leslie Kean

# Load libraries
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(janitor)
library(rstatix)
library(readr)
library(readxl)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/06_TCR_Diversity"))

# Clear environment variables
rm(list = ls())

# Load the tsv files for relapse cohort
rel_files <- list.files(
  paste0(my_wd, "/../../1801 data/Relapse TCRs/"),
  full.names = T
)
relapse_list <- lapply(rel_files, function(x) read_tsv(x, guess_max = 100000))
# Need to still look into the warnings and make sure they don't cause issues downstream
names(relapse_list) <- substr(basename(rel_files), 22, 24)

# Load the tsv files for remission cohort
rem_files <- list.files(
  paste0(my_wd, "/../../1801 data/Remission TCRs/"),
  full.names = T
)
remission_list <- lapply(rem_files, function(x) read_tsv(x, guess_max = 100000))
# Need to still look into the warnings and make sure they don't cause issues downstream
names(remission_list) <- substr(basename(rem_files), 20, 22)

# Concatenate 2 lists
concatenated_list <- c(remission_list, relapse_list)

# Remove the outliers (less than 10k or more than 100k rows)
rows_per_sample <- unlist(lapply(concatenated_list, nrow))
# Visualize
plot(sort(rows_per_sample), log = "y", ylab = "Number of rows")
abline(h = c(10000, 100000), col = "red")
# Subset
within_range <- between(rows_per_sample, 10000, 100000)
subset_list <- concatenated_list[within_range]
# Check that the minimum is >10k and the maximum is <100k
summary(unlist(lapply(subset_list, nrow)))

#Check for only Tac-MTX group
#Easiest way to do it
# Add the clinical info table
rel_tb <- read_xlsx(paste0(
  my_wd,
  "/../../1801 data/Clinical tables/relapse cohort tcr data.xlsx"
))
rem_tb <- read_xlsx(paste0(
  my_wd,
  "/../../1801 data/Clinical tables/remission cohort tcr data.xlsx"
))

# Merge 2 tables
clinical_tb <- rbind(rel_tb, rem_tb)
tac <- subset(clinical_tb, treatment_arm == "Tac/MTX")
# Merge calculations with clinical merged table
t1 <- merge(clinical_tb, merged)
#select only TAC-MTX group
tac <- subset(clinical_tb, treatment_arm == "Tac/MTX")
tac <- tac$Patient_ID
formatted_numbers <- ifelse(
  tac >= 10 & tac < 100,
  sprintf("%03d", tac),
  as.character(tac)
)

subset_list1 <- subset_list[names(subset_list) %in% formatted_numbers]


# Peter's addition 240702 --------------------------------------------------------------------------
# Example expansion of a single dataframe
# df <- subset_list[[1]]
# # Visualize
# df$templates %>% sort %>% plot(log = "y")
# plot(sort(df$templates), log = "y")
# # Expand
# expanded_df <- df %>% uncount(templates)
# expanded_df <- uncount(df, templates)
# # Check
# sum(df$templates) == nrow(expanded_df)

# pdf("6.5_Clone_sizes.pdf")
# lapply(subset_list, function(x) plot(sort(x$templates), log = "y", ylim = c(1, 42000)))
# dev.off()
# End of Peter's addition -------------------------------------------------------------------------

# Use lapply to expand all dataframes in the list
expanded_list <- lapply(subset_list1, function(x) uncount(x, templates))

# Check the expansion visually by plotting the original vs. the expanded number of rows for both lists
rows_per_sample2 <- unlist(lapply(expanded_list, nrow))
rows_per_sample_subset <- rows_per_sample[names(rows_per_sample2)]
plot(
  x = rows_per_sample_subset,
  y = rows_per_sample2,
  log = "xy",
  xlab = "Original row number",
  ylab = "After expanding"
)

# Use Ksenia's function from 6.1 to calculate indices
.diversityCall <- function(data) {
  shannon <- .shannon(data[, "Freq"])
  inv_simpson <- .invsimpson(data[, "Freq"])
  norm_entropy <- .normentropy(data[, "Freq"])
  gini_simpson <- .ginisimpson(data[, "Freq"])
  chao1 <- .chao1(data[, "Freq"])
  ACE <- .ACE(data[, "Freq"])
  out <- c(shannon, inv_simpson, norm_entropy, gini_simpson, chao1, ACE)
  return(out)
}

.shannon <- function(p) {
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(-sum(p * log(p)))
}
.normentropy <- function(p) {
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(-sum(p * log(p)) / log(length(p)))
}
.invsimpson <- function(p) {
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(1 / sum(p^2))
}
.ginisimpson <- function(p) {
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(1 - sum(p^2))
}
.chao1 <- function(p) {
  n1 <- sum(p == 1)
  n2 <- sum(p == 2)
  S_obs <- length(p)
  # Chao1 index calculation
  if (n1 > 1 && n2 > 0) {
    chao1 <- S_obs + (n1 * (n1 - 1)) / (2 * (n2 + 1))
  } else {
    # In cases where n1 <= 1 or n2 == 0, Chao1 is undefined
    chao1 <- NA
  }
  return(chao1)
}
.ACE <- function(p) {
  q <- 10
  S_abund <- sum(p > q)
  rare_data <- p[p <= q]
  S_rare <- length(rare_data)
  n_rare <- sum(rare_data)

  # Calculate C_ACE
  C_ACE <- sum(p) / n_rare

  # Calculate gamma
  gamma <- 0
  for (i in seq_len(q)) {
    f_i <- sum(rare_data == i)
    gamma <- gamma + (1 - i / q)^f_i
  }

  # Calculate ACE
  ACE <- S_abund + (S_rare / C_ACE) + (1 - C_ACE) * gamma
  return(ACE)
}

# Define the function
compute_diversity = function(df_list, cloneCall, n.boots) {
  # df_list <- combined.sc
  # cloneCall <- "CTstrict"
  # n.boots <- 1000
  mat = NULL
  min_n = 10000 # we are doing 500 here instead of minumum cell number across the samples ### UNDER CONSTRUCTION
  for (i in seq_along(df_list)) {
    data = df_list[[i]]
    mat_a = NULL
    for (j in seq(seq_len(n.boots))) {
      x = slice_sample(data, n = min_n)
      y = as.data.frame(table(x[, cloneCall]))
      sample = .diversityCall(y)
      mat_a = rbind(mat_a, sample)
    }
    mat_a[is.na(mat_a)] = 0
    mat_b = colMeans(mat_a)
    mat_b = as.data.frame(t(mat_b))
    mat = rbind(mat, mat_b)
  }
  colnames(mat) = c(
    "shannon",
    "inv.simpson",
    "norm.entropy",
    "gini.simpson",
    "chao1",
    "ACE"
  )
  mat$Patient_ID = names(df_list)

  return(mat)
}

# Calculate the diversity with the function, keep in mind only works with lists
m <- compute_diversity(
  expanded_list,
  cloneCall = "rearrangement",
  n.boots = 1000
)
View(m)
#m$inv.simpson <- round(m$inv.simpson, 2)
m_selected <- m[, c("Patient_ID", "inv.simpson")]

# Also convert Simpson clonality from Adaptive results to inverse Simpson index
extracted_columns <- lapply(expanded_list, function(df) {
  df[1, "sample_simpson_clonality"]
})
t_df <- data.frame(
  "Patient_ID" = names(extracted_columns),
  "simpson_clonality" = unlist(extracted_columns),
  row.names = NULL
)
t_df$inverse_simpson_index <- 1 / (t_df$simpson_clonality^2)
#t_df$inverse_simpson_index <- round(t_df$inverse_simpson_index, 2)
#t_df$simpson_clonality <- round(t_df$simpson_clonality, 2)

# Export the table
write.csv(t_df, "calculated_inv.simpson.csv")

# Compare 2 values by merging into same table (calculated vs converted)
merged <- inner_join(m_selected, t_df, by = "Patient_ID") %>%
  arrange(Patient_ID)

# Remove leading 0's in the patient IDs to merge with the clinical dataset
merged <- merged %>% mutate(Patient_ID = gsub("^0+", "", Patient_ID))

# Add the clinical info table
rel_tb <- read_xlsx(paste0(
  my_wd,
  "/../../1801 data/Clinical tables/relapse cohort tcr data.xlsx"
))
rem_tb <- read_xlsx(paste0(
  my_wd,
  "/../../1801 data/Clinical tables/remission cohort tcr data.xlsx"
))

# Merge 2 tables
clinical_tb <- rbind(rel_tb, rem_tb)

# Merge calculations with clinical merged table
t1 <- merge(clinical_tb, merged)

# Save the data
write.csv(t1, "ultimate_table.csv")

# Make a categorical variant according to the relapse point
# t1 <-  t1 %>%
#   mutate(Relapse_type = case_when(Relapse.Months.Post.Transplant <6 ~ "early_relapse",
#                                   Relapse.Months.Post.Transplant >=6 ~ "late_relapse"))

# Remove people that saying death and survival less then 11
t1_filtered <- t1 %>%
  filter(
    !(cohort == "remission" &
      Overall.Survival.Months.Post.Transplant < 11 &
      Overall.Survival.Post.Transplant.Event == "Death")
  )

# Also remove the patients in relapse cohort who already relapsed before 100 day (118,171,222,426)
t1_filtered_final <- t1_filtered %>%
  filter(!(cohort == "relapse" & Relapse.Days.Post.Transplant < 100))

# Save the data
write.csv(t1_filtered_final, "filtered_table.csv")

# Compare the Simpson indexes that was calculated differently
p <- ggplot(data = merged, aes(x = inverse_simpson_index, y = inv.simpson)) +
  geom_point(aes(color = Patient_ID), size = 4) +
  theme_bw() +
  ylab("Inverse Simpson Index/calculated") +
  xlab("Inverse Simpson Index/converted") +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 15,
      color = "black"
    ),
    axis.title.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15, color = "black"),
    legend.key.size = unit(3, "mm"),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
p

# compare the simpson indexes, logaritmic
p1 <- ggplot(
  data = t1_filtered_final,
  aes(x = inverse_simpson_index, y = inv.simpson)
) +
  geom_point(size = 4) +
  theme_bw() +
  ylab("Inverse Simpson Index/calculated") +
  xlab("Inverse Simpson Index/converted") +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 15,
      color = "black"
    ),
    axis.title.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15, color = "black"),
    legend.key.size = unit(3, "mm"),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_x_log10()
p1

# Visualize using different parameters from the clinical table
p2 <- t1_filtered_final %>%
  filter(!duplicated(Patient_ID)) %>%
  ggplot(aes(x = cohort, y = inv.simpson)) +
  geom_boxplot() +
  geom_jitter(aes(color = treatment_arm), size = 4) +
  scale_color_manual(values = c("skyblue1", "salmon")) +
  theme_bw() +
  ylab("Inverse Simpson Index") +
  theme(strip.text = element_text(size = 14, color = "black", face = "bold")) +
  theme(
    aspect.ratio = 1.5,
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 15,
      color = "black"
    ),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15, color = "black"),
    legend.key.size = unit(3, "mm"),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  stat_compare_means(
    aes(group = cohort),
    method = "wilcox",
    method.args = list(var.equal = T),
    label = "p.format",
    label.x = 1.5,
    label.y = 300,
    tip.length = 1,
    size = 4
  )
p2

# Visualize for each treatment arm
p3 <- t1_filtered_final %>%
  #filter(Indicator.for.GRFS.event..1.0.==0) %>%
  ggplot(aes(x = Event.for.Primary.Endpoint..GRFS, y = inv.simpson)) +
  geom_boxplot() +
  geom_jitter(aes(color = cohort), size = 4) +
  scale_color_manual(values = c("skyblue1", "salmon")) +
  theme_bw() +
  ylab("Inverse Simpson Index") +
  theme(strip.text = element_text(size = 14, color = "black", face = "bold")) +
  theme(
    aspect.ratio = 1.5,
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 15,
      color = "black"
    ),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15, color = "black"),
    legend.key.size = unit(3, "mm"),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
+stat_compare_means(
  aes(group = cohort),
  method = "wilcox",
  method.args = list(var.equal = T),
  label = "p.format",
  label.x = 1.5,
  label.y = 300,
  tip.length = 1,
  size = 4
)
p3

pdf("1801results_type.pdf", width = 8, height = 8)
p2
dev.off()
