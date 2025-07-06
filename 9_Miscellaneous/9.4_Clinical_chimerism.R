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
p1 <- ggplot(df2, aes(x = Cohort, y = Chimerism, color = Cohort)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.5) +
  coord_cartesian(ylim = c(0, 101)) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    size = 0.4,
    color = "#00000080"
  ) +
  stat_compare_means(
    aes(group = Cohort),
    method = "wilcox.test",
    label.y = 50,
    label.x = 1.25,
    size = 3,
    label = "p.format",
    show.legend = F
  ) +
  labs(y = "Donor chimerism (clinical)") +
  scale_color_manual(
    values = c("long-term-remission" = "#546fb5", "relapse" = "#e54c35")
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 2,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank()
  )

# View
p1

# Save it as PDF
pdf("9.4_Clinical_chimerism.pdf", width = 4, height = 4)
p1
dev.off()
