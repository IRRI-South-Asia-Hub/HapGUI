pca_plot <- function(pcafile){
  
  #PCA PLOT
  m <- read.csv(pcafile, header = TRUE)
  colnames(m) <- m[1,]
  m <- m[-c(1),]
  colnames(m)[1] <- "ID"
  
  pop <- read.delim("Edited all population.txt")
  
  names(pop) <- c("IRIS.ID","Name","SUBPOPULATION","COUNTRY","IRGC.NO")
  comb_p <- merge(m, pop, by.x = "ID", by.y = "IRIS.ID", all.x = T)
  
  
  comb_p[,2:5] <- lapply(comb_p[,2:5],as.numeric)
  comb_p[,2:5] <- lapply(comb_p[,2:5], round, 2)
  
  theme_cus <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                    strip.background=element_blank(),
                    axis.text.y=element_text(colour="black", face = "bold"),
                    axis.text.x=element_text(colour="black", face = "bold",angle = 90),
                    axis.ticks=element_line(colour="black"),
                    strip.text.x = element_text(colour="black", face = "bold"),
                    title = element_text(color = "black", size = 12, family = "Arial", face = "bold"),
                    plot.margin=unit(c(1,1,1,1),"line"),
                    legend.title = element_text(color = "black", size = 12, family = "Arial", face = "bold"),
                    legend.text = element_text(colour = "black", size = 12,family = "Arial"),
                    legend.key = element_rect(fill = "white", color = NA))
  p1 <- ggplot(comb_p) + geom_point(aes(x=PC1,y=PC2, color=SUBPOPULATION)) +
    theme_cus + theme(legend.position = "none")+
    scale_color_manual(breaks = c("temp", "aus", "trop", "indx", "ind1A", "ind1B", "admix", "ind2", "subtrop", "japx", "ind3", "aro"),
                       values = c("#F1948A", "#D7BDE2", "#E6B0AA", "#D4E6F1", "#2E86C1", "#2E86C1", "#5499C7", "#CCD1D1", "#85C1E9", "#CD6155", "#E74C3C", "#A9CCE3", "#A9DFBF"))
  
  p2 <- ggplot(comb_p) + geom_point(aes(x=PC1,y=PC3, color=SUBPOPULATION)) +
    theme_cus + theme(legend.position = "none")+
    scale_color_manual(breaks = c("temp", "aus", "trop", "indx", "ind1A", "ind1B", "admix", "ind2", "subtrop", "japx", "ind3", "aro"),
                       values = c("#F1948A", "#D7BDE2", "#E6B0AA", "#D4E6F1", "#2E86C1", "#2E86C1", "#5499C7", "#CCD1D1", "#85C1E9", "#CD6155", "#E74C3C", "#A9CCE3", "#A9DFBF"))
  
  p3 <- ggplot(comb_p) + geom_point(aes(x=PC2,y=PC3, color=SUBPOPULATION)) +
    theme_cus + theme(legend.position = "none")+
    scale_color_manual(breaks = c("temp", "aus", "trop", "indx", "ind1A", "ind1B", "admix", "ind2", "subtrop", "japx", "ind3", "aro"),
                       values = c("#F1948A", "#D7BDE2", "#E6B0AA", "#D4E6F1", "#2E86C1", "#2E86C1", "#5499C7", "#CCD1D1", "#85C1E9", "#CD6155", "#E74C3C", "#A9CCE3", "#A9DFBF"))
  
  ar <- ggarrange(p1,ggarrange(p2,p3,labels = c("B", "C"),ncol = 2,nrow = 1),
                  ncol = 1,nrow = 2,labels = "A",
                  common.legend = T, legend = "bottom")
  ar
  
  ggsave("pca_plot_2D.png", ar)
}
