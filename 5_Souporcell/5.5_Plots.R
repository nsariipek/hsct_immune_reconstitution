# Nurefsan Sariipek, 231103
# Generating plots using souporcell information
# Updated at 240418

# Load the libraries
library(ggalluvial)
library(tidyverse)
library(janitor)
library(ggforce)
library(RColorBrewer)
library(randomcoloR)
library(ggrepel)
library(fossil)
library(ggpubr)
library(readxl)

# Empty environment
rm(list=ls())

# Set the working directory 
# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Souporcell/"

# For Peter:
# my_wd <- "~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Souporcell/"
# Load the saved dataframe that contains the information from the previous part
# For Nurefsan
# For cohorts 1-2
combined_df <- read_csv(paste0(my_wd,file = "/results/cohort1-2_souporcell.csv"))
# For cohort 3
combined_df <- read_csv(paste0(my_wd,file="/results/cohort3_souporcell.csv"))

# Subset and summarise as you like to visualize 

######################### Plot 1 ##########################
# Plot 1, visualize souporcell results per patient
# Select the patient you want to visualize 
pt1 <- subset(combined_df, subset = patient_identity=="pt01")
pt2 <- subset(combined_df, subset = patient_identity=="pt02")
pt5 <- subset(combined_df, subset = patient_identity=="pt05")
pt6 <- subset(combined_df, subset = patient_identity=="pt06")
pt8 <- subset(combined_df, subset = patient_identity=="pt08")
pt9 <- subset(combined_df, subset = patient_identity=="pt09")
pt10 <- subset(combined_df, subset = patient_identity=="pt10")
pt11 <- subset(combined_df, subset = patient_identity=="pt11")
pt12 <- subset(combined_df, subset = patient_identity=="pt12")


# Re-order timepoints for better visualization
ordered_timepoints= c("pre_transplant", "remission","relapse")
pt1$status= factor(pt1$status, levels = ordered_timepoints)

# Calculate the proportion for donor/host cells
pt1 <- pt1 %>%
  mutate(assignment = gsub("host", "recipient",assignment))%>%
 group_by(status) %>% reframe(tabyl(assignment)) %>%
  mutate(percent = percent*100) 

# Plot

plot1 <- ggplot(data = pt1,
                aes(x = status,y = percent,fill= assignment))+  
                geom_bar(stat = "identity", position = "stack")+
               theme_bw() +
  scale_fill_manual(values = c("donor"= "tomato1", "recipient"="royalblue1"),name = "Assignment")+
  labs(y = "Percent of Cells") +
  theme(aspect.ratio = 0.75,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,face="plain", size=18, color="black"), 
        axis.text.y = element_text(face="plain", size=16, color="black"),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        plot.title = element_text(size=20, face="plain"),
        strip.text = element_text(size=12, face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=16)) +
        theme(panel.grid = element_blank()
              )

plot1
# Save as a pdf
pdf("pt1.pdf", width = 6, height = 8 )
plot1
dev.off()


# Select only T and NK cells
combined_df_T <- subset(combined_df, subset = celltype %in% c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Bright NK cells","CD56 Dim NK cells")) 

# Select only 6mo remission cells and MNC libraries only to prevent any skewing
combined_df_T_rem6mo <- subset(combined_df_T, subset = id %in% c("P01.1Rem","P01.2Rem","P02.1Rem","P04.1Rem", "P05.1Rem","P06.1Rem","P07.1Rem","P08.1Rem"))
 
##################### Plot 2 #############################
# Plot 2, Post-transplant 3-6 months remission samples alluvial visuliazing

# Summarize the dataframe
meta_summary <- 
  combined_df_T_rem6mo  %>% group_by(celltype,cohort,id,assignment) %>%
  mutate(cohort = gsub("Remission cohort", "Non-relapsed", gsub("Relapse cohort", "Relapsed", cohort)), cohort = factor(cohort, levels = c("Non-relapsed", "Relapsed")))  %>%
  summarize(n_cells = n())

# Add the color palette that Peter has used in his UMAP 
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pdf("230724_Colorwheel.pdf")
pie(rep(1,74), col=col_vector)
dev.off()

mycol_tib <- tribble(~celltype, ~color,
                     "CD4 Naïve", col_vector[23],
                     "CD4 Memory", col_vector[28],
                     "Treg", col_vector[11],
                     "CD8 Naïve", col_vector[67],
                     "CD8 Memory", col_vector[63],
                     "CD8 Effector", col_vector[58],
                     "CD8 Terminally Exhausted", col_vector[60],
                     "γδ T lymphocytes", col_vector[1],
                     "NK T cells", col_vector[46],
                     "CD56 Bright NK cells", col_vector[61],
                     "CD56 Dim NK cells", col_vector[52])
mycol <- mycol_tib$color
names(mycol) <- mycol_tib$celltype

# Visualize

plot2 <- ggplot(data = meta_summary,
       aes(axis1 = cohort , axis2 = assignment, axis3 = celltype,
           y = n_cells)) +
  scale_x_discrete(limits = c("Cohort","Origin","Celltype"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = celltype)) + #change fill to change the flow
  geom_stratum(aes(fill=celltype), color="black") + 
  geom_text(stat = "stratum", aes(label =paste(after_stat(stratum))), size =4) +
  scale_fill_manual(values = mycol) + 
  theme_bw()+
  ylab("Number of cells")+
  theme(aspect.ratio = 0.5,
    axis.text.x = element_text(face="bold", size=16, color="black"), 
        axis.title.x = element_text(face="bold", size=16, color="black"),
    axis.title.y = element_text(face="bold", size=16, color="black"),
        axis.text.y = element_text(face="plain", size=12, color="black"),
        strip.text = element_text(size=10, face="bold"),
        plot.title = element_text(size=20, face="plain"),
    legend.position = "none") +
  ggtitle("Post-transplant 3-6 months remission samples")

plot2
# Save as a pdf
pdf("3-6moalluvial.pdf", width = 18, height = 9)
plot2
dev.off()

######################### Plot 3 ##########################
# Plot 3, see how accurate souporcell results comparing with chimerism information

# Load the dataframe contains comparision between chimerism and souporcell
df <- read_excel(paste0(my_wd,"AnalysisNurefsan/Souporcell/outputs/similarity.xlsx"))
View(df)

ids = c("P01.1Rem","P01.1RemT", "P01.2Rem" , "P03.1Rem" , "P04.1Rem",  "P04.1RemT", "P05.1Rem" , "P05.Rel" ,  "P05.RelT" ,"P06.1Rem",  "P07.1Rem",  "P08.1Rem" , "P08.2Rem" ,"P09.1Rem", "P09.Rel", "P09.RelT", "P10.1Rem", "P10.Rel", "P11.1Rem" , "P11.1RemT", "P12.1Rem" )
df$id = factor(df$id, levels= ids)

# Visualize
plot3 =
  df %>% 
  drop_na()%>%
  ggplot(aes(x = Souporcell, y = Chimerism, colour = id)) +
  geom_point(size=3.5, position = "jitter")+
  #scale_color_distiller(palette = "Spectral")+
  #geom_text(hjust=0, vjust=0)+ also add label= id inside of aes 
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  coord_cartesian(expand = FALSE, xlim = c(0, 1.1), ylim = c(0, 1.1))+
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.title = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"))+
  theme(legend.text=element_text(size=8),legend.title=element_text(size=10, hjust = 0.5))+
  ggtitle("Correlation of Souporcell Results \n with Clinical Chimerism Data")+
  theme(aspect.ratio=1)

plot3
# Save as a pdf
pdf("similarity.pdf", width = 5, height = 7.5)
plot3
dev.off()

