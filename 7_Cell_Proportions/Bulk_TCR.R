

library(ggpubr)

# For Nurefsan:
 my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"

# read excel file
 t <- read_csv(paste0(my_wd,"Bulk/IrepSeq.csv"))


 
 
   p1 <- t %>% mutate(Vial = factor(Vial, levels = c("1","2","3")))%>%
   ggplot(aes(x=Group, y=`Diversity Index`)) +
   geom_jitter(aes(color=Vial), size=4)+
     xlab("Separated population") +
   theme_pubr() +
   theme(strip.text = element_text(size = 14, color = "black", face="bold"),
         aspect.ratio = 1.5, 
         axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 15, color = "black"),
         axis.text.y = element_text(size = 15),
         axis.title.y = element_text(size = 15, color = "black"),
         legend.key.size = unit(3,"mm"),
         legend.position = "right",
         legend.text = element_text(size = 12)) +
   expand_limits(x = 0, y = 0) + 
   stat_summary(fun.data=mean_se, geom="errorbar", width=.5, linewidth=1) +
   stat_compare_means(aes(group = Group), method = "t.test", 
                      method.args = list(var.equal = T),
                      label = "p.format", label.x = 1, 
                      label.y= 4, tip.length = 1, size = 6, inherit.aes = TRUE) 
p1 

p1
pdf("BulkIrep.pdf", width = 5, height = 7.5)
p1
dev.off()
getwd()
