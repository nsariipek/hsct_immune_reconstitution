# Peter van Galen, 250206
# Save standard cell type colors for throughout the project

library(tidyverse)
library(ggsci) # for scale_color_igv

# Set working directory
setwd("~/TP53_ImmuneEscape/2_Annotate/")

# Delete environment variables
rm(list=ls())

# The igv panel has 51 colors
color_data <- data.frame(
  index = factor(1:51),
  color = rev(pal_igv("default")(51))
)

# Generate the bar plot
ggplot(color_data, aes(x = 1, y = index, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +  # Use the color values directly
  scale_y_discrete(labels = color_data$color) +
  theme_minimal() +
  labs(x = "Hex Color Codes", title = "IGV Color Palette") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

# Now assign to cell types
celltype_colors <- c(`Progenitors` = "#3B1B53FF",
  `Early Erythroids` = "#D60047FF",
  `Mid Erythroids` = "#924822FF",
  `Late Erythroids` = "#AE1F63FF",
  `Pro Monocytes` = "#99CC00FF",
  `Monocytes` = "#E4AF69FF",
  `Non Classical Monocytes` = "#7A65A5FF",
  `cDC` = "#5DB1DDFF",
  `pDC` = "#CDDEB7FF",
  `Pro B cells` = "#14FFB1FF",
  `Pre-B` = "#00991AFF",
  `B cells` = "#003399FF",
  `Plasma cells` = "#802268FF",
  `CD4 Naïve` = "#466983FF",
  `CD4 Effector Memory` = "#D58F5CFF",
  `CD4 Memory` = "#C75127FF",
  `Treg` = "#FFC20AFF",
  `CD8 Naïve` = "#33CC00FF",
  `CD8 Effector` = "#612A79FF",
  `CD8 Memory` = "#0099CCFF",
  `CD8 Exhausted` = "#CE3D32FF",
  `γδ T` = "#D595A7FF",
  `NK T` = "#5050FFFF",
  `Adaptive NK` = "#1A0099FF",
  `CD56 Bright NK` = "#00D68FFF",
  `CD56 Dim NK` = "#008099FF",
  `Cycling T-NK cells` = "#F0E685FF",
  `UD1` = "#A9A9A9FF",
  `UD2` = "#837B8DFF",
  `UD3` = "#5A655EFF",
  `pal_igv_unused01` = "#BA6338FF",
  `pal_igv_unused02` = "#CC9900FF",
  `pal_igv_unused03` = "#99CC00FF",
  `pal_igv_unused04` = "#E7C76FFF",
  `pal_igv_unused05` = "#CC9900FF",
  `pal_igv_unused06` = "#00CC99FF",
  `pal_igv_unused07` = "#4775FFFF",
  `pal_igv_unused08` = "#00CC33FF",
  `pal_igv_unused09` = "#0A47FFFF",
  `pal_igv_unused10` = "#990033FF",
  `pal_igv_unused11` = "#991A00FF",
  `pal_igv_unused12` = "#996600FF",
  `pal_igv_unused13` = "#809900FF",
  `pal_igv_unused14` = "#749B58FF",
  `pal_igv_unused15` = "#339900FF",
  `pal_igv_unused16` = "#009966FF",
  `pal_igv_unused17` = "#660099FF",
  `pal_igv_unused18` = "#990080FF",
  `pal_igv_unused19` = "#6BD76BFF",
  `pal_igv_unused20` = "#FFD147FF",
  `pal_igv_unused21` = "#FF1463FF")

# Check
all(celltype_colors %in% color_data$color)
all(color_data$color %in% celltype_colors)

# Save
write.table(data.frame(celltype = names(celltype_colors), color = celltype_colors, stringsAsFactors = FALSE), "../celltype_colors.txt", sep = "\t", row.names = FALSE, quote = FALSE)