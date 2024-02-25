
system("softwares_external/plink2 --bfile  --pca 5 --nonfounders --out pc")

#PLINK
pca_file <- paste0(infam, "-pc.eigenvec")
val_file <- paste0(infam, "-pc.eigenval")
m <- read.table(pca_file, header = F)
cat(blue(paste("Loading",pca_file,"...\n")))
names(m)[1:5] <- c("FID","IID","PC1","PC2","PC3")

phe <- read.delim(pheno, header = T)
cat(blue(paste("Loading",pheno,"...\n")))
names(phe) <- c("Name","IRIS.ID","SUBPOPULATION","COUNTRY")

comb_p <- merge(m, phe, by.x = "FID", by.y = "IRIS.ID", all.x = T)

p1 <- ggplot(comb_p) + geom_point(aes(x=PC1,y=PC2, color=SUBPOPULATION)) +
  theme_cus + theme(legend.position = "none")+
  scale_color_manual(breaks = c("indx","ind2","ind1B","ind1A","ind3",
                                "aus","japx","temp",
                                "trop","subtrop","admix","aro" ),
                     values = c( "#56B4E9", "#56B4E9","#0072B2", "#293352",
                                 "#00AFBB", "#000000", "#E69F00","#F0E442",
                                 "#D55E00", "#CC79A7","#999999", "#52854C"))
p1 <- p1 + theme_Publication()

ggsave(paste0(trait,"_pca.jpeg"), p1, width = 8, height = 8)

p2 <- ggplot(comb_p) + geom_point(aes(x=PC1,y=PC3, color=SUBPOPULATION)) +
  theme_cus + theme(legend.position = "none")+
  scale_color_manual(breaks = c("indx","ind2","ind1B","ind1A","ind3",
                                "aus","japx","temp",
                                "trop","subtrop","admix","aro" ),
                     values = c( "#56B4E9", "#56B4E9","#0072B2", "#293352",
                                 "#00AFBB", "#000000", "#E69F00","#F0E442",
                                 "#D55E00", "#CC79A7","#999999", "#52854C"))
p3 <- ggplot(comb_p) + geom_point(aes(x=PC2,y=PC3, color=SUBPOPULATION)) +
  theme_cus + theme(legend.position = "none")+
  scale_color_manual(breaks = c("indx","ind2","ind1B","ind1A","ind3",
                                "aus","japx","temp",
                                "trop","subtrop","admix","aro" ),
                     values = c( "#56B4E9", "#56B4E9","#0072B2", "#293352","#00AFBB",
                                 "#000000", "#E69F00","#F0E442",
                                 "#D55E00", "#CC79A7","#999999", "#52854C"))

ar <- ggarrange(p1,ggarrange(p2,p3,labels = c("B", "C"),ncol = 2,nrow = 1),
                ncol = 1,nrow = 2,labels = "A", 
                common.legend = T, legend = "bottom")
print(ar)
cat(green(paste0("Saving PCA plot for first 4 components in ",pca_plot, "...\n")))
ggsave(pca_plot, ar, width = 8, height = 6)

value <- read.table(val_file, header = F)
var <- (value/sum(value))*100
var$V2 <- c(1:nrow(value))

pov <- ggplot(var) + geom_line(aes(x = V2, y = V1), color="red") +
  geom_point(aes(x = V2, y = V1),shape = 8) + 
  theme_cus + xlab("Principal Components") + ylab("PoV")

cat(green(paste0("Saving PCA PoV of first 4 components in pca_pov.pdf", 
                 pov_plot,"...\n")))
ggsave(pov_plot, pov, width = 8, height = 6)