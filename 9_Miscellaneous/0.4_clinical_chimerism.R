# Load packages
library(tidyverse)
library(ggpubr)

# Import the data
df <- read_csv("~/clinical_chimerism.csv") 
cohort_colors <- c("long-term-remission" = "#546fb5FF","relapse" = "#e54c35ff")

# Prepare the data
df <- df %>%
  mutate(
    Cohort = as.factor(Cohort),
    Chimerism = as.numeric(Chimerism)
  ) %>%
  filter(!is.na(Chimerism))

#Boxplot
ggplot(df, aes(x = Cohort, y = Chimerism, fill = Cohort)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.6, size = 1.5) +
  scale_fill_manual(values = cohort_colors) +
  ylim(0,100)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        ylab("Clinical Chimerism"),
        axis.title.x = element_blank())+
  stat_compare_means(method = "t.test", label = "p.format", label.y = 75,label.y.npc = "center", size = 4.5)
                     
# Save it as PDF
ggsave("chimerism_comparison.pdf", width = 6, height = 8)
