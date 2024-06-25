candidate_gene <- function(pos_infile,ld, genes_file){
  res <- read.csv(pos_infile)
  res$ID <- paste(res$Chromosome,"_",res$Marker.position..bp.)
  #res <- read.csv("PL_hap_dup.csv")
  ricegenes <- read.delim(genes_file)
  cangene_ind <-  res$ID
  candigene_ind <-list()
  ld <- as.numeric(ld)
  
  for (pos in cangene_ind) {

    chrom <- subset(res,res$ID%in% pos)
    qtn <- chrom$Marker.position..bp.
    star <- qtn-ld
    stop <- qtn+ld
    selgenes <- ricegenes[ricegenes$start >=star & ricegenes$start <=stop,]
    selgenes <- data.frame(subset(selgenes,selgenes$Chromosome %in% chrom$Chromosome))

    if(nrow(selgenes) > 0){
      selgenes$QTN <- c(pos)
      selgenes <- data.frame(selgenes)
      snp <- paste0("rs_",pos)
      candigene_ind[[snp]] <- selgenes
    }
  }    
  candidates <- bind_rows(candigene_ind)
  
  names(candidates)[1] <- "gene_id"
  meth <- unlist(strsplit(pos_infile,split = "_"))[1]
  write.csv(candidates,paste0(meth,"_candidate_genes.csv"),row.names = F)
  
  return(candidates)
}
