
# Pie Chart 
# Load libraries
library(tidyverse)

# Your survival colors
survival_colors <- c(
  "Non-relapsed" = "#4775FFFF",
  "Relapsed" = "#E64B35FF"
)

# Make lighter versions for WT
lighter_colors <- c("Relapsed_WT" = scales::alpha("#E64B35FF", 0.4),
  "Non-relapsed_WT" = scales::alpha("#4775FFFF", 0.4)
)

# Create the data
df <- tibble(
  tp53_status = c("TP53MT", "TP53WT", "TP53MT", "TP53WT"),
  outcome = c("Relapsed", "Relapsed", "Non-relapsed", "Non-relapsed"),
  n = c(8, 6, 6, 13)
) %>%
  mutate(
    label = paste(tp53_status, outcome, "\n(n=", n, ")", sep=""),
    fill_group = case_when(
      tp53_status == "TP53MT" ~ outcome,
      tp53_status == "TP53WT" ~ paste0(outcome, "_WT")
    )
  )

# Calculate positions
df <- df %>%
  arrange(outcome, tp53_status) %>%
  mutate(csum = cumsum(n),
         pos = csum - n/2)

# Plot
pie_chart <- ggplot(df, aes(x = "", y = n, fill = fill_group)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(y = pos, label = label), color = "black", size = 3) +
  scale_fill_manual(values = c(lighter_colors,survival_colors)) +
  theme_void() +
  theme(legend.position = "none")


pdf("Figure1_piechart.pdf", width = 8, height = 8)
pie_chart
dev.off()


  