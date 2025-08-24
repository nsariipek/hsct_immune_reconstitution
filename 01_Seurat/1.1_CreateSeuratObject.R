# Nurefsan Sariipek and Peter van Galen, updated on 250417
# Create Seurat object with metadata from CellRanger count matrices

# Load the needed libraries
library(tidyverse)
library(Seurat) # we are using version 5.1.0
library(googleCloudStorageR)
library(janitor)

# Set working directory (in Terra). Setting the home directory is an exception to my usual advice of using the folder with the script.
setwd("/home/rstudio/")
# Google VM:
#setwd("/home/unix/vangalen")

# Start with a clean slate
rm(list = ls())

# Parameters to interact with Google bucket, this part only needed for Terra
gcs_global_bucket("fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b")
# Check if you can list the objects. In Terra, you may need to authenticate using gcs_auth(). In VM, this did not work - hence the alternative function on line 65.
gcs_list_objects()

# Load matrices of all samples
Samples <- c(
  "P1013_MNC",
  "P1732_MNC",
  "P1953_MNC",
  "P2434_MNC",
  "P2599_CD3",
  "P2791_MNC",
  "P6174_CD3",
  "P9931_CD3",
  "P1195_MNC",
  "P1745_MNC",
  "P1964_MNC",
  "P2446_MNC",
  "P2599_MNC",
  "P2820_MIX",
  "P6174_MNC",
  "P9931_MNC",
  "P1285_MNC",
  "P1762_MIX",
  "P1972_CD3",
  "P2448_MNC",
  "P2621_CD3",
  "P2961_MNC",
  "P6244_CD3",
  "P1347_CD3",
  "P1764_MIX",
  "P1972_MNC",
  "P2517_MIX",
  "P2621_MNC",
  "P2977_MIX",
  "P6244_MNC",
  "P1347_MNC",
  "P1804_MNC",
  "P2220_MNC",
  "P2518_CD3",
  "P2645_MNC",
  "P2986_MNC",
  "P9185_CD3",
  "P1665_MIX",
  "P1811_CD3",
  "P2332_MNC",
  "P2518_MNC",
  "P2698_MIX",
  "P2988_MNC",
  "P9185_MNC",
  "P1671_MIX",
  "P1811_MNC",
  "P2379_CD3",
  "P25809_MNC",
  "P2737_CD3",
  "P3000_MIX",
  "P9355_MNC",
  "P1677_CD3",
  "P1817_MIX",
  "P2379_MNC",
  "P25802_CD3",
  "P2737_MNC",
  "P4618_MNC",
  "P9596_CD3",
  "P1677_MNC",
  "P1953_CD3",
  "P2408_MNC",
  "P25802_MNC",
  "P2745_MNC",
  "P5641_MNC",
  "P9596_MNC"
)

# If it already exists, make sure this folder is empty
dir.create("~/tmp")

# Function to download data and create Seurat object (for Terra)
process_sample <- function(Sample) {
  #Sample <- Samples[1]
  print(Sample)

  # Create path
  sample_path <- paste0(Sample, "/sample_filtered_feature_bc_matrix/")

  # Download files
  files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  for (file in files) {
    gcs_get_object(
      object_name = paste0(sample_path, file),
      saveToDisk = file.path("/home/rstudio/tmp", file)
    )
  }

  # Read data and create Seurat object
  data <- Read10X(data.dir = "/home/rstudio/tmp")
  seu <- CreateSeuratObject(counts = data, project = Sample)

  # Remove downloaded files
  unlink(file.path("/home/rstudio/tmp", files))

  # Create Seurat object
  return(seu)
}

# Function to download data and create Seurat object (for VM)
process_sample <- function(Sample) {
  #Sample <- Samples[1]
  print(Sample)

  # Create path
  sample_path <- paste0(Sample, "/sample_filtered_feature_bc_matrix/")

  # Download files
  files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  for (file in files) {
    system(paste0(
      "gsutil cp gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/",
      sample_path,
      file,
      " /home/unix/vangalen/tmp"
    ))
  }

  # Read data and create Seurat object
  data <- Read10X(data.dir = "/home/unix/vangalen/tmp")
  seu <- CreateSeuratObject(counts = data, project = Sample)

  # Remove downloaded files
  unlink(file.path("/home/unix/vangalen/tmp", files))

  # Create Seurat object
  return(seu)
}

# Use lapply to populate the list with Seurat objects
seu_ls <- lapply(Samples, process_sample)

# Create the Seurat object combining each sample, add.cell.ids prevents duplicate cell identifiers
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)], add.cell.ids = Samples)
# The following takes a while and requires ~100 GB memory
seu <- JoinLayers(seu)

# Add QC Metrics
seu <- PercentageFeatureSet(seu, "^MT-", col.name = "percent_mito")
seu <- PercentageFeatureSet(seu, "^RB[SL]", col.name = "percent_ribo")

# Check metadata
head(seu@meta.data)

# Visualize QC metrics as a violin plot. First make the object 10x smaller to speed this up.
seu_subset_for_visualization <- subset(
  seu,
  cells = colnames(seu)[seq(1, ncol(seu), by = 10)]
)

# Inspect plots, grouped by patient id, for any potential outlier samples
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")
VlnPlot(
  seu_subset_for_visualization,
  group.by = "orig.ident",
  features = feats,
  pt.size = 0.1,
  ncol = 2,
  alpha = 0.3
) +
  NoLegend()
ggsave(
  "~/hsct_immune_reconstitution/1_Seurat/1.1_FeatureViolins.pdf",
  width = 25,
  height = 10
)

FeatureScatter(
  seu_subset_for_visualization,
  "nCount_RNA",
  "nFeature_RNA",
  group.by = "orig.ident",
  pt.size = 0.5,
  shuffle = T,
  plot.cor = F
) +
  theme(aspect.ratio = 1) +
  geom_vline(xintercept = 250, col = "black") +
  geom_hline(yintercept = 500, col = "black")
ggsave(
  "~/hsct_immune_reconstitution/1_Seurat/1.1_FeatureScatter.pdf",
  width = 10,
  height = 6
)

gc()

# Count the cells before filtering
cell_numbers <- as.data.frame(
  seu$orig.ident %>% tabyl %>% select(".", "n") %>% rename(n_prefilter = n)
)

# Filter data by QC thresholds based on the plots above
seu <- subset(
  seu,
  subset = nFeature_RNA > 250 & nCount_RNA > 500 & percent_mito < 20
)

# Check how many cells were filtered
t1 <- as.data.frame(
  seu$orig.ident %>% tabyl %>% select(".", "n") %>% rename(n_postfilter = n)
)
cellnumbers_final <- merge(t1, cell_numbers, by = ".") %>%
  mutate(percent = n_postfilter / n_prefilter)
colnames(cellnumbers_final)[1] <- "orig.id"
write.csv(
  cellnumbers_final,
  "~/hsct_immune_reconstitution/1_Seurat/cellnumbers_final.csv"
)

# Add variables to metadata
# Add library type
seu$library_type <- case_when(
  grepl("CD3", seu$orig.ident) ~ "CD3",
  grepl("MNC", seu$orig.ident) ~ "MNC",
  grepl("MIX", seu$orig.ident) ~ "MIX"
)

# Add cohort (whether the patient relapsed or not)
seu$cohort <- case_when(
  grepl(
    "2446|25802|2645|1972|2220|2621|9185|2599|1665|1745|1817|2408|2988|1762|2698|2791|2977|2986|1671|2517|2820|2961|3000",
    seu$orig.ident
  ) ~
    "long-term-remission",
  grepl(
    "9596|2737|2379|2434|2518|4618|6174|9931|1953|25809|1677|5641|1732|1811|1195|1347|1285|6244|9355|1013|1764|1804|1964|2332|2448|2745",
    seu$orig.ident
  ) ~
    "relapse"
)

# Add cohort details
seu$cohort_detail <- case_when(
  grepl("9596|2737|2379|2434|2518|4618|6174|9931|1953|25809", seu$orig.ident) ~
    "1-Relapse",
  grepl("2446|25802|2645|1972|2220|2621|9185|2599", seu$orig.ident) ~
    "1-Long-term-remission",
  grepl("1677|5641|1732|1811|1195|1347|1285|6244|9355|1013", seu$orig.ident) ~
    "1-Early-relapse",
  grepl("1764|1804|1964|2332|2448|2745", seu$orig.ident) ~ "2-Relapse",
  grepl(
    "1665|1745|1817|2408|2988|1762|2698|2791|2977|2986|1671|2517|2820|2961|3000",
    seu$orig.ident
  ) ~
    "2-Long-term-remission"
)

# Add sample status
seu$sample_status <- case_when(
  grepl("9596|2379|4618|9355|2446|1972|1677|1195|5641", seu$orig.ident) ~
    "pre-transplant",
  grepl(
    "25809|2434|2518|6174|9931|1013|25802|2645|2220|2621|9185|2599|1732|1285|6244|1764|1804|2332|2745|1665|1745|1817|2408|2988|1762|2698|2791|2977|2986|1671|2517|2820|2961|3000",
    seu$orig.ident
  ) ~
    "remission",
  grepl("1964|2448|2737|1953|1811|1347", seu$orig.ident) ~ "relapse"
)

# Add patient number
seu$patient_id <- case_when(
  grepl("2446|25802|2645", seu$orig.ident) ~ "P01", # previously P01
  grepl("1972|2220|2621", seu$orig.ident) ~ "P02", # previously P02
  grepl("9185", seu$orig.ident) ~ "P03", # previously P03
  grepl("2599", seu$orig.ident) ~ "P04", # previously P04
  grepl("9596|25809|2737", seu$orig.ident) ~ "P20", # previously P05
  grepl("2379|2434", seu$orig.ident) ~ "P21", # previously P06
  grepl("2518", seu$orig.ident) ~ "P22", # previously P07
  grepl("4618|6174|9931|1953", seu$orig.ident) ~ "P23", # previously P08
  grepl("1677|1732|1811", seu$orig.ident) ~ "P30", # previously P09
  grepl("1195|1285|1347", seu$orig.ident) ~ "P31", # previously P10
  grepl("5641|6244", seu$orig.ident) ~ "P32", # previously P11
  grepl("9355|1013", seu$orig.ident) ~ "P33", # previously P12
  grepl("1665", seu$orig.ident) ~ "P07", # previously P13
  grepl("1745", seu$orig.ident) ~ "P05", # previously P14
  grepl("1817", seu$orig.ident) ~ "P08", # previously P15
  grepl("2408", seu$orig.ident) ~ "P09", # previously P16
  grepl("2988", seu$orig.ident) ~ "P06", # previously P17
  grepl("1762", seu$orig.ident) ~ "P10", # previously P18
  grepl("2698", seu$orig.ident) ~ "P11", # previously P19
  grepl("2791", seu$orig.ident) ~ "P12", # previously P20
  grepl("2977", seu$orig.ident) ~ "P13", # previously P21
  grepl("2986", seu$orig.ident) ~ "P14", # previously P22
  grepl("1671", seu$orig.ident) ~ "P15", # previously P23
  grepl("2517", seu$orig.ident) ~ "P16", # previously P24
  grepl("2820", seu$orig.ident) ~ "P17", # previously P25
  grepl("2961", seu$orig.ident) ~ "P18", # previously P26
  grepl("3000", seu$orig.ident) ~ "P19", # previously P27
  grepl("1764", seu$orig.ident) ~ "P24", # previously P28
  grepl("1804", seu$orig.ident) ~ "P25", # previously P29
  grepl("1964", seu$orig.ident) ~ "P28", # previously P30
  grepl("2332", seu$orig.ident) ~ "P27", # previously P31
  grepl("2448", seu$orig.ident) ~ "P29", # previously P32
  grepl("2745", seu$orig.ident) ~ "P26"
) # previously P33

# Add a unique sample identifier to combine MNC and CD3 libraries
seu$sample_id <- case_when(
  grepl("2446", seu$orig.ident) ~ "P01_Pre", # previously P01
  grepl("25802", seu$orig.ident) ~ "P01_Rem1", # previously P01
  grepl("2645", seu$orig.ident) ~ "P01_Rem2", # previously P01
  grepl("1972", seu$orig.ident) ~ "P02_Pre", # previously P02
  grepl("2220", seu$orig.ident) ~ "P02_Rem1", # previously P02
  grepl("2621", seu$orig.ident) ~ "P02_Rem2", # previously P02
  grepl("9185", seu$orig.ident) ~ "P03_Rem", # previously P03
  grepl("2599", seu$orig.ident) ~ "P04_Rem", # previously P04
  grepl("9596", seu$orig.ident) ~ "P20_Pre", # previously P05
  grepl("25809", seu$orig.ident) ~ "P20_Rem", # previously P05
  grepl("2737", seu$orig.ident) ~ "P20_Rel", # previously P05
  grepl("2379", seu$orig.ident) ~ "P21_Pre", # previously P06
  grepl("2434", seu$orig.ident) ~ "P21_Rem", # previously P06
  grepl("2518", seu$orig.ident) ~ "P22_Rem", # previously P07
  grepl("4618", seu$orig.ident) ~ "P23_Pre", # previously P08
  grepl("6174", seu$orig.ident) ~ "P23_Rem1", # previously P08
  grepl("9931", seu$orig.ident) ~ "P23_Rem2", # previously P08
  grepl("1953", seu$orig.ident) ~ "P23_Rel", # previously P08
  grepl("1677", seu$orig.ident) ~ "P30_Pre", # previously P09
  grepl("1732", seu$orig.ident) ~ "P30_Rem", # previously P09
  grepl("1811", seu$orig.ident) ~ "P30_Rel", # previously P09
  grepl("1195", seu$orig.ident) ~ "P31_Pre", # previously P10
  grepl("1285", seu$orig.ident) ~ "P31_Rem", # previously P10
  grepl("1347", seu$orig.ident) ~ "P31_Rel", # previously P10
  grepl("5641", seu$orig.ident) ~ "P32_Pre", # previously P11
  grepl("6244", seu$orig.ident) ~ "P32_Rem", # previously P11
  grepl("9355", seu$orig.ident) ~ "P33_Pre", # previously P12
  grepl("1013", seu$orig.ident) ~ "P33_Rem", # previously P12
  grepl("1665", seu$orig.ident) ~ "P07_Rem", # previously P13
  grepl("1745", seu$orig.ident) ~ "P05_Rem", # previously P14
  grepl("1817", seu$orig.ident) ~ "P08_Rem", # previously P15
  grepl("2408", seu$orig.ident) ~ "P09_Rem", # previously P16
  grepl("2988", seu$orig.ident) ~ "P06_Rem", # previously P17
  grepl("1762", seu$orig.ident) ~ "P10_Rem", # previously P18
  grepl("2698", seu$orig.ident) ~ "P11_Rem", # previously P19
  grepl("2791", seu$orig.ident) ~ "P12_Rem", # previously P20
  grepl("2977", seu$orig.ident) ~ "P13_Rem", # previously P21
  grepl("2986", seu$orig.ident) ~ "P14_Rem", # previously P22
  grepl("1671", seu$orig.ident) ~ "P15_Rem", # previously P23
  grepl("2517", seu$orig.ident) ~ "P16_Rem", # previously P24
  grepl("2820", seu$orig.ident) ~ "P17_Rem", # previously P25
  grepl("2961", seu$orig.ident) ~ "P18_Rem", # previously P26
  grepl("3000", seu$orig.ident) ~ "P19_Rem", # previously P27
  grepl("1764", seu$orig.ident) ~ "P24_Rem", # previously P28
  grepl("1804", seu$orig.ident) ~ "P25_Rem", # previously P29
  grepl("1964", seu$orig.ident) ~ "P28_Rel", # previously P30
  grepl("2332", seu$orig.ident) ~ "P27_Rem", # previously P31
  grepl("2448", seu$orig.ident) ~ "P29_Rel", # previously P32
  grepl("2745", seu$orig.ident) ~ "P26_Rem"
) # previously P33

# Add TP53 mutation status
seu$TP53_status <- case_when(
  grepl(
    "P01|P02|P03|P04|P05|P06|P20|P21|P22|P23|P30|P31|P32|P33",
    seu$patient_id
  ) ~
    "MUT",
  grepl(
    "P07|P08|P09|P10|P11|P12|P13|P14|P15|P16|P17|P18|P19|P24|P25|P26|P27|P28|P29",
    seu$patient_id
  ) ~
    "WT"
)

# Add the time point to show the sample time as months after tx
seu$timepoint <- case_when(
  grepl("2446", seu$orig.ident) ~ 0,
  grepl("25802", seu$orig.ident) ~ 3,
  grepl("2645", seu$orig.ident) ~ 6,
  grepl("1972", seu$orig.ident) ~ 0,
  grepl("2220", seu$orig.ident) ~ 5,
  grepl("2621", seu$orig.ident) ~ 24,
  grepl("9185", seu$orig.ident) ~ 12,
  grepl("2599", seu$orig.ident) ~ 3,
  grepl("9596", seu$orig.ident) ~ 0,
  grepl("25809", seu$orig.ident) ~ 3,
  grepl("2737", seu$orig.ident) ~ 9,
  grepl("2379", seu$orig.ident) ~ 0,
  grepl("2434", seu$orig.ident) ~ 3,
  grepl("2518", seu$orig.ident) ~ 3,
  grepl("4618", seu$orig.ident) ~ 0,
  grepl("6174", seu$orig.ident) ~ 3,
  grepl("9931", seu$orig.ident) ~ 12,
  grepl("1953", seu$orig.ident) ~ 24,
  grepl("1677", seu$orig.ident) ~ 0,
  grepl("1732", seu$orig.ident) ~ 1,
  grepl("1811", seu$orig.ident) ~ 3,
  grepl("1195", seu$orig.ident) ~ 0,
  grepl("1285", seu$orig.ident) ~ 1.5,
  grepl("1347", seu$orig.ident) ~ 3,
  grepl("5641", seu$orig.ident) ~ 0,
  grepl("6244", seu$orig.ident) ~ 1.5,
  grepl("9355", seu$orig.ident) ~ 0,
  grepl("1013", seu$orig.ident) ~ 1.5,
  grepl("1665", seu$orig.ident) ~ 3,
  grepl("1745", seu$orig.ident) ~ 3,
  grepl("1817", seu$orig.ident) ~ 3,
  grepl("2408", seu$orig.ident) ~ 3,
  grepl("2988", seu$orig.ident) ~ 3,
  grepl("1762", seu$orig.ident) ~ 3,
  grepl("2698", seu$orig.ident) ~ 3,
  grepl("2791", seu$orig.ident) ~ 3,
  grepl("2977", seu$orig.ident) ~ 3,
  grepl("2986", seu$orig.ident) ~ 3,
  grepl("1671", seu$orig.ident) ~ 3,
  grepl("2517", seu$orig.ident) ~ 3,
  grepl("2820", seu$orig.ident) ~ 3,
  grepl("2961", seu$orig.ident) ~ 3,
  grepl("3000", seu$orig.ident) ~ 3,
  grepl("1764", seu$orig.ident) ~ 3,
  grepl("1804", seu$orig.ident) ~ 3,
  grepl("1964", seu$orig.ident) ~ 3,
  grepl("2332", seu$orig.ident) ~ 3,
  grepl("2448", seu$orig.ident) ~ 3,
  grepl("2745", seu$orig.ident) ~ 3
)

# Convert variables to factors
seu$orig.ident <- as.factor(seu@meta.data$orig.ident)
seu$library_type <- factor(
  seu@meta.data$library_type,
  levels = c("MNC", "CD3", "MIX")
)
seu$cohort <- factor(
  seu@meta.data$cohort,
  levels = c("long-term-remission", "relapse")
)
seu$cohort_detail <- factor(
  seu@meta.data$cohort_detail,
  levels = c(
    "1-Long-term-remission",
    "1-Relapse",
    "1-Early-relapse",
    "2-Long-term-remission",
    "2-Relapse"
  )
)
seu$sample_status <- factor(
  seu@meta.data$sample_status,
  levels = c("pre-transplant", "remission", "relapse")
)
seu$patient_id <- as.factor(seu@meta.data$patient_id)
seu$sample_id <- as.factor(seu@meta.data$sample_id)
# Check that there are no missing values
sapply(seu@meta.data, function(x) sum(is.na(x)))

# Save (this takes a while, you can monitor progress (growing file size) in the Terminal - it's about 2.2 Gb in the end)
saveRDS(
  seu,
  file = "~/hsct_immune_reconstitution/AuxiliaryFiles/250417_MergedSeuratObject.rds"
)
