qc_linux <- function(choice_geno,phefile,ip_dir, theme_cus) {
  
  system(paste0(ip_dir, "/plink2 --vcf marker.vcf --maf 0.05 --export vcf --out marker"))

  system(paste0(ip_dir, "/tassel5/run_pipeline.pl -Xmx12g -fork1 -vcf marker.vcf -homozygous -export marker -exportType Hapmap -runfork1"))
  
  system(paste0(ip_dir, "/plink2 --vcf marker.vcf --freq --out marker"))
  system("bcftools annotate --set-id '%CHROM\\_%POS' -o marker_updated.vcf marker.vcf")
  system(paste0(ip_dir, "/plink2 --vcf marker.vcf --pca --read-freq marker.afreq --out pca"))
  
  
  mar <- read.delim2(paste0(dir,"/marker.hmp.txt"), header = FALSE)
  mar <- mar[, -c(2, 5:11)]
  pca1 <- read.delim2(paste0(dir,"/pca.eigenvec"), header = FALSE)
  
  if (choice_geno == "option1") {
    names <- mar[1, 4:ncol(mar)]
    names2 <- sapply(names, function(x) unlist(strsplit(x, "_IRIS"))[1])
    
    pca1[c(2:nrow(pca1)),1] <- names2
  } else{
    names2 <- mar[1, 4:ncol(mar)]
  }
  
  colnames(mar) <- c("rs#", "chrom", "pos", names2)
  mar <- mar[-1,]

  write.csv(mar, paste0(dir,"/marker.csv"), row.names = FALSE)

  
  colnames(pca1) <- c("ID", paste0("PC", 1:(ncol(pca1) - 1)))
  pca1 <- pca1[-c(1),]
  pca1[,2:11] <- lapply(pca1[,2:11], as.numeric)

  write.csv(pca1, file = paste0(dir,"/pca.csv"), row.names = FALSE)
  
  # comb_p <- merge(pca1, phefile, by = "ID")
  comb_p <- merge(pca1, phefile, by = "ID")
  comb_p$subpop <- factor(comb_p$subpop)
  
  # p1 <- ggplot(comb_p, aes(x = PC1, y = PC2, color = subpop)) +
  #   geom_jitter()+theme_cus

  # p1 <- ggplot(comb_p, aes(x = PC1, y = PC2, color = subpopulation)) + geom_jitter()+theme_cus
  # p1 <- ggplot(comb_p, aes(x = PC1, y = PC2, color = subpop)) + geom_jitter() + theme_cus
  
  if (choice_geno == "option1") {
    p1 <- ggplot(comb_p, aes(x = PC1, y = PC2, color = subpop)) + 
      geom_jitter() + theme_cus
  } else if (choice_geno == "option2") {
    p1 <- ggplot(comb_p, aes(x = PC1, y = PC2, color = subpopulation)) + 
      geom_jitter() + theme_cus
  }
  print(colnames(comb_p))
  
  ggsave(paste0(dir,"/pca_plot_2D.png"), p1, width = 12, height = 10, dpi = 600)
  
  return(mar)
}
