candidate_gene <- function(pos_infile,ld, gff_file){


  res <- read.csv(pos_infile)
  colnames(res)[1:2] <- c("Chromosome", "Position")
  print(head(res))   
  res$ID <- paste(res$Chromosome, "_", res$Position, sep = "")
  
  #add code for gff file
  gff <- read.table(gff_file, sep = "\t", quote = "", stringsAsFactors = FALSE)
  genes <- gff[gff$V3 == "gene", ]

  # g1 <- genes %>% separate(V9, into = c("key", "value"), sep = ";")
  # g1 <- g1[,-ncol(g1)]
  # g1 <- g1 %>% separate(key, into = c("key", "value"), sep = ":")
  
  g1 <- genes %>%
    separate(V9, into = c("key", "value"), sep = ";", extra = "merge", fill = "right") %>%
    separate(key, into = c("key", "value"), sep = ":", extra = "merge", fill = "right")
  
  g1 <- g1[,!names(g1) %in% c("key","V2","V3","V6", "V8")]
  colnames(g1) <- c("Chromosome", "start", "stop", "strand", "gene_id")


  g1 <- g1 %>% distinct(gene_id, .keep_all = TRUE)

  ## Roman to int
  roman_vals <- c(I=1, II=2, III=3, IV=4, V=5, VI=6, VII=7, VIII=8, IX=9, X=10)
  g1$Chromosome <- as.character(g1$Chromosome)
  g1$Chromosome <- ifelse(toupper(g1$Chromosome) %in% names(roman_vals),
                          as.character(roman_vals[toupper(g1$Chromosome)]),
                          g1$Chromosome)
  
  g11 <- g1[grepl("^[0-9]+$", g1$Chromosome),]
  g11$Chromosome <- as.character(g11$Chromosome)
  res$Chromosome <- as.character(res$Chromosome)

  rgene <- g11
  print(head(rgene)) 
  if (nrow(rgene) == 0) {
    stop("Error: No genes found in GFF file. Please check the format.")
  }

  candidates_list <- list()
  ld <- as.numeric(ld)

  for (i in 1:nrow(res)) {
    chrom <- res$Chromosome[i]
    pos <- as.numeric(res$Position[i])
    qtn_id <- res$ID[i]
    
    gene_hits <- g1 %>%
      filter(Chromosome == chrom) %>%
      filter((start >= (pos - ld) & start <= (pos + ld)) |
               (stop >= (pos - ld) & stop <= (pos + ld)) |
               (start <= (pos - ld) & stop >= (pos + ld))) %>%
      mutate(QTN = qtn_id)
    
    if (nrow(gene_hits) > 0) {
      candidates_list[[qtn_id]] <- gene_hits
    }
  }
  candidates <- bind_rows(candidates_list)
  write.csv(candidates, "candidate_genes.csv", row.names = FALSE)
  
  locus <- candidates %>%
    dplyr::select(gene_id) %>%
    dplyr::distinct() %>%
    dplyr::rename(df2 = gene_id)
  write.csv(locus, "locus.csv", row.names = FALSE)

  
  return(candidates)
}
