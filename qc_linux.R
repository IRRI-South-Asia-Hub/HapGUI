qc_linux <- function(phefile, dir){

  phe <- read.csv(phefile, header = T)
  colnames(phe) <- c("Accessions","trait")
  phe$Accessions <- gsub(pattern = "IRIS ",replacement = "IRIS_",
                         phe$Accessions)
  
  out1 <- phe[,c(1,1)]
  write.table(out1, "id.txt", col.names = F, row.names = F,
              quote = F)

  geno <- "server2"
  plink <- "plink2"
  tassel <- "run_pipeline.pl"
  
  print("all the things:")
  print(geno)
  print(plink)
  print(tassel)
  
  system(command = paste0(plink," --bfile ",geno," --keep id.txt --export vcf --out marker"))

  system(command = paste0(tassel," -fork1 -vcf marker.vcf -export marker -exportType Hapmap"))

  mar <- read.delim2(file = file.path(dir,"marker.hmp.txt"),header = F)

  mar <- mar[,-c(2,5:11)]

  names <- mar[1,4:ncol(mar)] #Extracting IRIS id's from row one
  #names <- as.data.frame(sub("IRIS_", "IRIS-", names, fixed = TRUE))
  names2 <- data.frame(matrix(vector(), nrow = 1, ncol = ncol(names)))

  for (i in c(1:ncol(names))) {
    names2[1,i] <- unlist(strsplit(names[1,i],"_IRIS"))[1]
    names2[1,i] <- sub("IRIS_","IRIS-",names2[1,i])
  }

  #data processing
  col <- c("rs#","chrom","pos") #making new vector with these colnames

  #replace dots with -

  colnames(mar) <- c(col,names2) #combining the vector in the frame

  mar <- mar[-c(1),]

  write.csv(mar, file = file.path(dir,"marker.csv"),row.names = FALSE) #write processed data to CSV

  # PCA ---------------------------------------------------------------------

  system(command = paste0(plink," --vcf ",dir,"/marker.vcf --freq --out ",dir,"/freq"))

  system(command = paste0(plink, " --vcf ",dir,"/marker.vcf --pca --read-freq ",
                          dir,"/freq.afreq --out ",dir,"/pca"))

  ##reading PCA data
  pca1 <- read.delim2(file = file.path(dir,"pca.eigenvec"),header = FALSE) #1

  names <- as.data.frame(sub("IRIS_", "IRIS-", pca1$V1, fixed = TRUE))

  for (i in c(2:nrow(pca1))) {
    pca1[i,1] <- unlist(strsplit(names[i,1],"_"))[1]
  }

  colnames(pca1) <- "<PCA>"
  colnames(pca1)[2:ncol(pca1)] <- ""
  pca1[1,1] <- "<ID>"

  #writing PCA data to CSV
  write.csv(pca1, file = file.path(dir,"pca.csv"), row.names = FALSE)
  
  #PCA PLOT
  m <- read.csv(file.path(dir,"pca.csv"), header = TRUE)
  colnames(m) <- m[1,]
  m <- m[-c(1),]
  colnames(m)[1] <- "ID"
  
  file_path <- system.file("","Edited all population.txt", package = "HaploGUI")
  pop <- read.delim(file_path)
  
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
  
  ggsave(file.path(dir,"pca_plot_2D.png"), ar)
}
