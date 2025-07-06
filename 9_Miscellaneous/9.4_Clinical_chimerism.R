# Nurefsan Sariipek and Peter van Galen, 250705
# Plot clinical chimerism during remission

# Load packages
library(tidyverse)
library(ggpubr)

# Set working directory
setwd(
  "~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Miscellaneous"
)

# Clear environment variables
rm(list = ls())

# Import the data
df <- read_csv("9.4_Clinical_chimerism.csv")

# Prepare the data
df2 <- df %>%
  mutate(
    Cohort = as.factor(Cohort),
    Chimerism = as.numeric(Chimerism)
  )

# Boxplot
ggplot(df2, aes(x = Cohort, y = Chimerism, fill = Cohort)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.6, size = 1.5) +
  scale_fill_manual(
    values = c("long-term-remission" = "#546fb5", "relapse" = "#e54c35")
  ) +
  ylim(0, 101) +
  labs(y = "Clinical chimerism") +
  theme_bw() +
  theme(
    aspect.ratio = 2,
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",
    label.x = 1.4,
    label.y = 75,
    size = 4.5
  )

# Save it as PDF
ggsave("9.4_Clinical_chimerism.pdf", width = 5, height = 4)
