# Nurefsan Sariipek and Peter van Galen, 250705
# Make the swimmer plot using the sample information for Figure 1

library(tidyverse)
library(readxl)

# Set working directory
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Miscellaneous")

# Clear environment variables
rm(list=ls())

# Load the table has sample timepoint info
df <- read_excel("9.2_Timepoints.xlsx")

# Wrangle the data
df2 <- df %>% select(-Note) %>%
  mutate(
    Patient_id = factor(Patient_id, levels = rev(unique(Patient_id))),
    Sample_type = case_when(
      str_detect(Sample_id, "pre") ~ "Pre",
      str_detect(Sample_id, "tx") ~ "Transplant",
      str_detect(Sample_id, "rem") ~ "Remission",
      str_detect(Sample_id, "rel") ~ "Relapse",
      TRUE ~ "Other"
    ),
    Patient_numeric = as.numeric(Patient_id))

# Extend lines for patients whose last sample is a Remission
max_days <- 24 * 30.44

extended_df <- df2 %>%
  group_by(Patient_id) %>%
  filter(`Timepoint(days)` == max(`Timepoint(days)`)) %>%
  filter(Sample_type == "Remission") %>%
  mutate(`Timepoint(days)` = max_days)

df_line <- bind_rows(df2, extended_df)

# Create background rectangles
rect_df <- df2 %>%
  group_by(Patient_id) %>%
  arrange(`Timepoint(days)`, .by_group = TRUE) %>%
  summarise(
    xmin = min(`Timepoint(days)`, na.rm = TRUE),
    last_tp = max(`Timepoint(days)`, na.rm = TRUE),
    last_type = Sample_type[which.max(`Timepoint(days)`)],
    .groups = "drop"
  ) %>%
  mutate(
    xmax = ifelse(last_type == "Remission", max_days, last_tp),
    Patient_numeric = as.numeric(factor(Patient_id, levels = rev(unique(df$Patient_id)))),
    ymin = Patient_numeric - 0.3,
    ymax = Patient_numeric + 0.3
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
  geom_rect(data = rect_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray85", alpha = 0.5) +
  geom_rect(
    aes(xmin = 90, xmax = 180, ymin = -Inf, ymax = Inf),
    fill ="#E30B5C", alpha = 0.2
  ) +
  annotate("text",
           x = 135,
           y = max(df2$Patient_numeric, na.rm = TRUE) + 1,
           label = "3–6 months",
           size = 4,
           fontface = "italic",
           color = "gray30") +
  geom_line(data = df_line,
            aes(x = `Timepoint(days)`, y = Patient_numeric, group = Patient_id),
            color = "gray30", linewidth = 0.8) +
  geom_point(data = df2,
             aes(x = `Timepoint(days)`, y = Patient_numeric, color = Sample_type),
             size = 3) +
  scale_y_continuous(name = "Patient ID",
    breaks = df2$Patient_numeric,
    labels = df2$Patient_id,
    expand = c(0.01, 0.01)
  ) +
  scale_x_continuous(
    name = "Months After Transplant",
    breaks = seq(0, 24 * 30.44, by = 3 * 30.44),
    labels = function(x) as.character(round(x / 30.44)),
    limits = c(min(df2$`Timepoint(days)`, na.rm = TRUE) - 10, 24 * 30.44),
    expand = c(0.01, 0.01)
  ) +
  scale_color_manual(values = sample_type_colors) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
  aspect.ratio = 1.5)

pdf("9.2_Swimmer_plot.pdf", width = 8,height = 10)
swimmer
dev.off()
