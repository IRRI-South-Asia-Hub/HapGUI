qc_windows <- function(phefile, ip_dir, choose_geno){
  
  phe <- read.csv(phefile, header = T)
  colnames(phe) <- c("Accessions","trait")
  phe$Accessions <- gsub(pattern = "IRIS ",replacement = "IRIS_",
                         phe$Accessions)
  
  out1 <- phe[,c(1,1)]
  write.table(out1, file = "id.txt", col.names = F, 
              row.names = F,quote = F)
  
  if (choose_geno == "Default"){
    genofile = "marker_main"
  }
  
  if (choose_geno == "External"){
    genofile = "extra_marker"
  }
  
  system(command = paste0(ip_dir, "/plink2 --bfile ", genofile, " --keep id.txt --export vcf --out marker_temp"))
  
  system(command = paste0(ip_dir, "/plink2 --vcf marker_temp.vcf --maf 0.05 --export vcf --out marker"))
  
  system(command = paste0(ip_dir,"/tassel-5-standalone/run_pipeline.bat -fork1 -vcf marker.vcf -homozygous -export marker -exportType Hapmap"))
  
  mar <- read.delim2(file = "marker.hmp.txt",header = F)
  
  mar <- mar[,-c(2,5:11)]
  
  names <- mar[1,4:ncol(mar)] #Extracting IRIS id's from row one
  #names <- as.data.frame(sub("IRIS_", "IRIS-", names, fixed = TRUE))
  names2 <- data.frame(matrix(vector(), nrow = 1, ncol = ncol(names)))
  
  for (i in c(1:ncol(names))) {
    names2[1,i] <- unlist(strsplit(names[1,i],"_IRIS"))[1]
    #names2[1,i] <- sub("IRIS_","IRIS-",names2[1,i])
  }
  
  #data processing
  col <- c("rs#","chrom","pos") #making new vector with these colnames
  
  #replace dots with -
  
  colnames(mar) <- c(col,names2) #combining the vector in the frame
  
  mar <- mar[-c(1),]
  
  write.csv(mar, file = "marker.csv",row.names = FALSE) #write processed data to CSV
  
  #Et-GWAS
  system(paste0(ip_dir,"/plink2 --bfile ", genofile ," --keep id.txt --make-bed --out marker"))
  
  # PCA ---------------------------------------------------------------------
  
  system(command = paste0(ip_dir, "/plink2 --vcf marker.vcf --freq --out marker"))
  
  system(command = paste0(ip_dir, "/plink2 --vcf marker.vcf --pca --read-freq marker.afreq --out pca"))
  
  ##reading PCA data
  pca1 <- read.delim2(file = "pca.eigenvec",header = FALSE) #1
  
  names <- as.data.frame(pca1$V1)
  
  for (i in c(2:nrow(pca1))) {
    pca1[i,1] <- unlist(strsplit(names[i,1],"_IRIS"))[1]
  }
  
  colnames(pca1) <- "<PCA>"
  colnames(pca1)[2:ncol(pca1)] <- ""
  pca1[1,1] <- "<ID>"
  
  #writing PCA data to CSV
  write.csv(pca1, file = "pca.csv", row.names = FALSE)
  
  #PCA PLOT
  m <- pca1
  colnames(m) <- m[1,]
  m <- m[-c(1),]
  colnames(m)[1] <- "ID"
  
  pop <- read.delim("Edited all population.txt")
  
  names(pop) <- c("IRIS.ID","Name","Subpopulation","COUNTRY","IRGC.NO")
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
  p1 <- ggplot(comb_p,aes(x=PC1,y=PC2, color=Subpopulation)) + 
    geom_jitter() +
    theme_cus +
    scale_color_manual(breaks = c("temp", "aus", "trop", "indx", "ind1A", "ind1B", "admix", "ind2", "subtrop", "japx", "ind3", "aro"),
                       values = c("#F1948A", "#D7BDE2", "#E6B0AA", "#D4E6F1", "#2E86C1", "#2E86C1", "#5499C7", "#CCD1D1", "#85C1E9", "#CD6155", "#E74C3C", "#A9CCE3", "#A9DFBF"))
  
  p2 <- ggplot(comb_p,aes(x=PC1,y=PC3, color=Subpopulation)) +
    geom_jitter()+
    theme_cus + theme(legend.position = "none")+
    scale_color_manual(breaks = c("temp", "aus", "trop", "indx", "ind1A", "ind1B", "admix", "ind2", "subtrop", "japx", "ind3", "aro"),
                       values = c("#F1948A", "#D7BDE2", "#E6B0AA", "#D4E6F1", "#2E86C1", "#2E86C1", "#5499C7", "#CCD1D1", "#85C1E9", "#CD6155", "#E74C3C", "#A9CCE3", "#A9DFBF"))
  
  p3 <- ggplot(comb_p,aes(x=PC2,y=PC3, color=Subpopulation)) +
    geom_jitter()+
    theme_cus + theme(legend.position = "none")+
    scale_color_manual(breaks = c("temp", "aus", "trop", "indx", "ind1A", "ind1B", "admix", "ind2", "subtrop", "japx", "ind3", "aro"),
                       values = c("#F1948A", "#D7BDE2", "#E6B0AA", "#D4E6F1", "#2E86C1", "#2E86C1", "#5499C7", "#CCD1D1", "#85C1E9", "#CD6155", "#E74C3C", "#A9CCE3", "#A9DFBF"))
  
  ar <- ggarrange(p1,ggarrange(p2,p3,labels = c("B", "C"),ncol = 2,nrow = 1),
                  ncol = 1,nrow = 2,labels = "A",
                  common.legend = T, legend = "bottom")
  ar
  
  ggsave("pca_plot_2D.png",plot = p1, width = 12, height = 10,
         dpi = 300, units = "cm")
  ggsave("pca_plot_all.png",plot = ar, width = 12, height = 10,
         dpi = 300, units = "cm")
  return(mar)
}
