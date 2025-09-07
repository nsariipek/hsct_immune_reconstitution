# Nurefsan Sariipek, Peter van Galen
# Figure 1D table

library(tidyverse)
library(readxl)

# Set working directory
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/10_Miscellaneous")

# Load table with sample mutation info
df <- read_excel("10.3_Mutation_table.xlsx")

# See number of mutations (for resutls section)
apply(select(df, -"Patient_id"), 2, function(x) sum(as.numeric(x)))

# Wrangle the df
df_long <- df %>%
  pivot_longer(cols = -Patient_id, names_to = "Gene", values_to = "Mutated")

df_long$Gene <- factor(
  df_long$Gene,
  levels = c(
    "TP53",
    "DNMT3A",
    "TET2",
    "ASXL1",
    "NPM1",
    "FLT3",
    "IDH1",
    "IDH2",
    "CK"
  )
)

p1 <- ggplot(
  df_long,
  aes(x = Gene, y = Patient_id, fill = as.factor(Mutated))
) +
  geom_tile(color = "black", size = 0.2) +
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

pdf("10.3_Mutation_plot.pdf", width = 2.7, height = 9.9)
p1
dev.off()
