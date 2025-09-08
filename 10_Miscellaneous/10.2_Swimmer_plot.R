# Nurefsan Sariipek and Peter van Galen, 250705
# Make the swimmer plot using the sample information for Figure 1

library(tidyverse)
library(readxl)
library(Seurat)

# Set working directory
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/10_Miscellaneous")

# Clear environment variables
rm(list = ls())

# Load the table has sample timepoint info
df <- read_excel("10.2_Timepoints.xlsx")

# Wrangle the data
df2 <- df %>%
  mutate(
    Timepoint_days = as.integer(gsub("Pre", -25, Timepoint_days)),
    y_placement = as.integer(factor(
      Patient_id,
      levels = rev(unique(Patient_id))
    ))
  )

# Check that the remission time point samples are consistent between Timepoints.xlsx and seu
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")
test_df <- cbind(
  swimmerfile = df2 %>%
    filter(between(Timepoint_days, 50, 200), Sample_type == "Remission") %>%
    pull(Sample_id) %>%
    sort,
  seuratfile = as_tibble(seu@meta.data) %>%
    filter(sample_status == "remission", timepoint %in% c(3, 5, 6)) %>%
    pull(sample_id) %>%
    as.character %>%
    unique %>%
    sort
)
identical(test_df[, "swimmerfile"], test_df[, "seuratfile"])

# Stats for manuscript text ---------------------------------------------------

# Long-term remission cohort follow up
df2 %>%
  filter(Cohort == "Long term remission") %>%
  group_by(Patient_id) %>%
  slice_max(order_by = Timepoint_days, n = 1) %>%
  mutate(Timepoint_months = Timepoint_days / 30.44) %>%
  pull(Timepoint_months) %>%
  (\(x) c(n = length(x), summary(x)))()
# --> "those who remained relapse-free during follow-up (n=19, median follow-up 36 months, range 8–83 months)"

# Time to relapse
df2 %>%
  filter(Cohort == "Relapse", Sample_type == "Relapse") %>%
  group_by(Patient_id) %>%
  filter(Timepoint_days == min(Timepoint_days)) %>%
  mutate(Timepoint_months = Timepoint_days / 30.44) %>%
  pull(Timepoint_months) %>%
  (\(x) c(n = length(x), summary(x)))()
# --> "those who experienced relapse (n=10, median time to relapse 5.5 months, range 2–39 months)"

# Remission sample collection time
df2 %>%
  filter(between(Timepoint_days, 50, 200), Sample_type == "Remission") %>%
  mutate(Timepoint_months = Timepoint_days / 30.44) %>%
  pull(Timepoint_months) %>%
  (\(x) c(n = length(x), summary(x)))()
# --> "We focused our analysis on remission samples collected ~3 months after transplant (n=27, median 3.3 months, range 2.2–5.9 months)"

# Time between sampling and relapse for patients with TP53 mutations
df2 %>%
  filter(
    Cohort == "Relapse",
    TP53_status == "MT",
    Sample_type %in% c("Remission", "Relapse")
  ) %>%
  select(Sample_id, Sample_type, Timepoint_days, Analyzed_cells)
sampling_to_relapse <- c(
  188 - 104,  #P20
  502 - 90,   #P21
  166 - 110,  #P22
  1201 - 116  #P23
) /
  30.44
summary(sampling_to_relapse)
# --> "The patients with TP53-mutated AML/MDS were in morphologic remission at this time and did not develop relapse until a median of 8.1 months after sampling (range 1.8–35.6 months)"

# Time between sampling and relapse for patients with >5% recipient HSPCs ~3 months after transplant
df2 %>%
  filter(
    Patient_id %in% c("P20", "P21", "P22", "P24", "P26", "P27"),
    Sample_type %in% c("Remission", "Relapse")
  ) %>%
  select(Sample_id, Sample_type, Timepoint_days, Analyzed_cells)
sampling_to_relapse <- c(
  188 - 104, #P20
  502 - 90,  #P21
  166 - 110, #P22
  168 - 88,  #P24
  133 - 104, #P26
  146 - 83   #P27
) /
  30.44
summary(sampling_to_relapse)
# --> "These six patients did not develop relapse until a median of 2.3 months after sampling (range 1.0–13.5 months). "

# Overall sample collection times (not used in the text)
df2 %>%
  filter(Timepoint_days != -25, Analyzed_cells == "Yes") %>%
  select(Patient_id, Timepoint_days) %>%
  print(n = 50)

###################

# Tibble for horizontal lines
df_line <- df2 %>%
  group_by(Patient_id) %>%
  filter(
    Timepoint_days == max(Timepoint_days) |
      Timepoint_days == min(Timepoint_days)
  ) %>%
  select(Patient_id, Timepoint_days, y_placement) %>%
  mutate(Timepoint_days = pmin(Timepoint_days, 63 * 30.44))

# Create background rectangles for remission (yellow) and relapse (red)
rem_rect <- df2 %>%
  group_by(Patient_id, y_placement) %>%
  summarise(
    xmin = 0,
    xmax = min(
      Timepoint_days[
        Sample_type %in% c("Last", "Relapse", "Death")
      ]
    )
  ) %>%
  mutate(
    xmax = pmin(xmax, 63 * 30.44),
    ymin = y_placement - 0.3,
    ymax = y_placement + 0.3
  )

rel_rect <- df2 %>%
  group_by(Patient_id, y_placement) %>%
  filter(Cohort != "Long term remission") %>%
  summarise(
    xmin = min(Timepoint_days[Sample_type == "Relapse"]),
    xmax = Timepoint_days[Sample_type == "Death"]
  ) %>%
  mutate(
    ymin = y_placement - 0.3,
    ymax = y_placement + 0.3
  )

# Define sample type colors
sample_type_colors <- c(
  "Pre" = "#A3BFD9",
  "Remission" = "#F6E06E",
  "Relapse" = "#8B0000",
  "Other" = "black"
)

# Create plot
swimmer <- ggplot() +
  geom_rect(
    aes(xmin = 68, xmax = 181, ymin = -Inf, ymax = Inf),
    fill = "#6d5689",
    alpha = 0.2
  ) +
  geom_line(
    data = df_line,
    aes(x = Timepoint_days, y = y_placement, group = Patient_id),
    color = "gray30",
    linewidth = 0.8
  ) +
  geom_rect(
    data = rem_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#ffd87e",
    alpha = 0.5
  ) +
  geom_rect(
    data = rel_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#f16858",
    alpha = 0.5
  ) +
  geom_point(
    data = filter(df2, Sample_type == "Death"),
    aes(x = Timepoint_days, y = y_placement),
    shape = 4,
    size = 3
  ) +
  geom_point(
    data = filter(df2, Timepoint_days != 0, Analyzed_cells == "Yes"),
    aes(x = Timepoint_days, y = y_placement, fill = Sample_type),
    stroke = 0.8,
    shape = 21,
    size = 3
  ) +
  scale_y_continuous(
    name = "Patient ID",
    breaks = unique(df2$y_placement),
    labels = unique(df2$Patient_id),
    expand = c(0.01, 0.01)
  ) +
  scale_x_continuous(
    name = "Months after transplant",
    breaks = seq(0, 63 * 30.44, by = 3 * 30.44),
    labels = function(x) as.character(round(x / 30.44)),
    limits = c(min(df2$Timepoint_days, na.rm = TRUE) - 10, 63 * 30.44),
    expand = c(0.01, 0.01)
  ) +
  scale_fill_manual(values = sample_type_colors) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    aspect.ratio = 0.5
  )

# View
swimmer

pdf("10.2_Swimmer_plot.pdf", width = 18, height = 8)
swimmer
dev.off()
