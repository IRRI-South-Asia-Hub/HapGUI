candidate_gene <- function(pos_infile,ld, gff_file){
  # pos_infile <- "/home/bandana/Documents/HapGUI_package/HapGUI-final_happhe/GWAS_results/1_Final result.csv"
  # gff_file <- "/media/bandana/DATA/HapGUI_revision/annotation_test/Oryza_sativa.IRGSP-1.0.60.gff3"
  res <- read.csv(pos_infile)
  res$ID <- paste(res$Chromosome,"_",res$Marker.position..bp.)
  #add code for gff file
  gff <- read.table(gff_file, sep="\t", quote="")
  genes <- gff[gff$V3 =="gene",]
  # g1 <- genes %>% separate(V9, into = c("key", "value"), sep = ";")
  # g1 <- g1[,-ncol(g1)]
  # g1 <- g1 %>% separate(key, into = c("key", "value"), sep = ":")

  g1 <- genes %>%
    separate(V9, into = c("key", "value"), sep = ";", extra = "merge", fill = "right") %>%
    separate(key, into = c("key", "value"), sep = ":", extra = "merge", fill = "right")

  g1 <- g1[,!names(g1) %in% c("key","V2","V3","V6", "V8")]
  colnames(g1) <- c("Chromosome", "start", "stop", "strand", "gene_id")
  g1 <- g1 %>% distinct(g1$gene_id, .keep_all = TRUE)
  g1 <- g1[,!names(g1) %in% c("g1$gene_id")]
  g11 <- g1[grepl("^[0-9]+$", g1$Chromosome),]
  g11$Chromosome <- as.character(g11$Chromosome)
  res$Chromosome <- as.character(res$Chromosome)

  rgene <- g11
  if (nrow(rgene) == 0) {
    stop("Error: No genes found in GFF file. Please check the format.")
  }

  cangene_ind <-  res$ID
  candigene_ind <-list()
  # ld <- as.numeric(50000)
  ld <- as.numeric(ld)

  for (pos in cangene_ind) {

    chrom <- subset(res,res$ID%in% pos)
    if (nrow(chrom) == 0) next
    qtn <- chrom$Marker.position..bp.
    star <- qtn-ld
    stop <- qtn+ld
    selgenes <- rgene[rgene$start >= star & rgene$start <= stop,]
    selgenes <- data.frame(subset(selgenes,selgenes$Chromosome %in% chrom$Chromosome))

    if(nrow(selgenes) > 0){
      # selgenes$QTN <- c(pos)
      selgenes$QTN <- pos
      selgenes <- data.frame(selgenes)
      snp <- paste0("rs_",pos)
      candigene_ind[[snp]] <- selgenes
    }
  }
  candidates <- bind_rows(candigene_ind)

  # names(candidates)[6] <- "gene_id"
  # meth <- unlist(strsplit(pos_infile,split = "_"))[1]
  write.csv(candidates,"candidate_genes.csv",row.names = F)

  return(candidates)
}
