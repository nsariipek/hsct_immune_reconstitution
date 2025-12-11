# Nurefsan Sariipek and Peter van Galen, 250528
# Save TCR data for all T cells

# Note: this script has to be run from Terra. Access to Google Cloud storage from other locations is finicky and it would be better to use `gsutil cp` instead, as in 1.1_CreateSeuratObject.R

# Load libraries
library(tidyverse)
library(Seurat)
library(googleCloudStorageR)
library(glue)
library(scRepertoire)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/06_TCR_Diversity"))

# Clear environment variables
rm(list = ls())

# Load Seurat object
seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")

# Favorite function
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Parameters to interact with Google bucket
gcs_global_bucket("fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b")
# Check if you can list the objects. Authenticate first.
gcs_auth()
gcs_list_objects()

# Define samples from Seurat object and compare to files in the bucket
Samples <- sort(as.character(unique(seu$orig.ident)))
vdj_folder_samples <- gcs_list_objects() %>%
  filter(grepl("vdj-t", name)) %>%
  pull(name) %>%
  cutf(., d = "/") %>%
  sort
identical(Samples, vdj_folder_samples) # should be TRUE

# Temporary directory to save downloaded files
tmp_dir <- "/home/rstudio/tmp"
dir.create(tmp_dir)

# Track successfully loaded samples
successful_samples <- c()

# Process each sample and assign to dynamically named variables
for (Sample in Samples) {
  #Sample <- Samples[1]
  print(Sample) # Log the sample being processed

  # Define file paths
  sample_path <- paste0(Sample, "/vdj-t/")
  file_name <- paste0(Sample, "_filtered_contig_annotations.csv")
  save_path <- file.path(tmp_dir, file_name)

  # Download the file and create variables
  tryCatch(
    {
      gcs_get_object(
        object_name = paste0(sample_path, "filtered_contig_annotations.csv"),
        saveToDisk = save_path
      )

      # Read the downloaded file into a variable named after the sample
      assign(Sample, read.csv(save_path), envir = globalenv())
      successful_samples <- c(successful_samples, Sample) # Track success
    },
    error = function(e) {
      message(glue("Error processing sample: {Sample}"))
    }
  )
}

# Check that they were all successful
Samples == successful_samples

# Construct the contig_list for successfully loaded samples
contig_list <- mget(successful_samples, envir = globalenv())

# Combine all the samples into a single object
TCR_combined_ls <- combineTCR(
  contig_list,
  samples = successful_samples
)

# Make it into a dataframe and select CTstrict
TCR_combined_df <- do.call(rbind, TCR_combined_ls)
TCR_combined_select_df <- select(TCR_combined_df, CTstrict)
rownames(TCR_combined_select_df) <- TCR_combined_df$barcode

# How many TCR calls have matching cells in the Seurat object? 99.2%
mean(rownames(TCR_combined_select_df) %in% colnames(seu))
TCR_combined_select_df <- TCR_combined_select_df[
  rownames(TCR_combined_select_df) %in% colnames(seu),
  ,
  drop = F
]

# Check: 37.1% of Seurat cells have a TCR
seu <- AddMetaData(seu, TCR_combined_select_df)
mean(!is.na(seu$CTstrict))

# Visualize the proportion of TCRs recovered per cell type --------------------

# Shortcut if you skipped the code above:
# seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Load colors
celltype_colors_df <- read.table(
  "../celltype_colors.txt",
  sep = "\t",
  header = T,
  stringsAsFactors = F,
  comment.char = ""
)
celltype_colors <- setNames(
  celltype_colors_df$color,
  celltype_colors_df$celltype
)

# Plot
as_tibble(seu@meta.data) %>%
  select(celltype, CTstrict) %>%
  group_by(celltype) %>%
  summarise(prop_TCR = mean(!is.na(CTstrict), na.rm = TRUE)) %>%
  ggplot(aes(x = celltype, y = prop_TCR, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  )

# Save plot
ggsave("6.1_TCR_calls.pdf", height = 4, width = 6)

# Save table
write_csv(
  as_tibble(TCR_combined_select_df, rownames = "cell"),
  file = "6.1_TCR_calls.csv.gz"
)


# In 2_Annotate-predict/2.2_Complete_Seurat_object.R, add TCR calls to Seurat object as follows:
tcr_calls <- read_csv("../06_TCR_Diversity/6.1_TCR_calls.csv.gz")
tcr_calls <- column_to_rownames(tcr_calls, var = "cell")
seu <- AddMetaData(seu, tcr_calls)
