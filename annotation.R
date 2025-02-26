vcf_annotation <- function(vcf_file, gff_file, ann_path, genome=NULL) {
  library(vcfR)
  library(reshape2)
  library(txdbmaker)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(BiocGenerics)
  library(S4Vectors)
  library(AnnotationDbi)
  library(VariantAnnotation)
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(rTASSEL)
  
  # Create output directory if it doesn't exist
  ann_path <- file.path(dir, "snpsift")
  dir.create(ann_path, recursive = TRUE, showWarnings = FALSE)
  
  # Read VCF and GFF files
  vcf_gr <- readVcf(vcf_file, genome = genome)
  gff <- read.table(gff_file, sep="\t", quote="")
  
  # Extract gene information from GFF
  genes <- gff[gff$V3 == "gene", ]
  g1 <- genes %>% separate(V9, into = c("key", "value"), sep = ";")
  g1 <- g1[,-ncol(g1)]
  g1 <- g1 %>% separate(key, into = c("key", "value"), sep = ":")
  g1 <- g1[,!names(g1) %in% c("key","V2","V3","V6", "V8")]
  colnames(g1) <- c("Chromosome", "start", "stop", "strand", "gene_id")
  g1 <- g1 %>% distinct(g1$gene_id, .keep_all = TRUE)
  g1 <- g1[,!names(g1) %in% c("g1$gene_id")]
  g11 <- g1[grepl("^[0-9]+(\\.[0-9]+)?$", g1$Chromosome), ]
  write.csv(g11, file.path(dir, "genes.csv"), row.names = FALSE)
  
  # Annotate variants
  txdb <- makeTxDbFromGFF(gff_file, format = "gff3")
  annotations <- locateVariants(vcf_gr, txdb, AllVariants())
  annotation_vcfgr <- as.data.frame(mcols(annotations))
  
  locations <- annotation_vcfgr$LOCATION
  gene_ids <- annotation_vcfgr$GENEID
  
  # Match annotations with VCF records
  variant_positions <- paste0(seqnames(annotations), ":", start(annotations))
  vcf_positions <- paste0(seqnames(vcf_gr), ":", start(vcf_gr))
  annotation_indices <- match(vcf_positions, variant_positions)
  
  locations <- locations[annotation_indices]
  gene_ids <- gene_ids[annotation_indices]
  
  # Add annotations to the VCF INFO field
  info(vcf_gr)$LOCATION <- locations
  info(vcf_gr)$GENEID <- gene_ids
  writeVcf(vcf_gr, file.path(dir, "annotated_variants.vcf"))
  
  annot_vcf <- read.vcfR(file.path(dir, "annotated_variants.vcf"))
  annot_info <- as.data.frame(annot_vcf@fix)
  geno_data <- as.data.frame(annot_vcf@gt)
  vcf_fin <- cbind(annot_info, geno_data)
  write.csv(vcf_fin, file.path(dir, "annotated_variants.csv"), row.names = FALSE)
  
  # Convert VCF to Hapmap format with MAF filtering
  tasObj <- readGenotypeTableFromPath(file.path(dir, "annotated_variants.vcf"))
  maf_vcfgr <- filterGenotypeTableSites(
    tasObj,
    siteMinCount = 0,
    siteMinAlleleFreq = 0.05
  )
  
  func.hmp <- read.delim("marker.hmp.txt", header = TRUE)
  
  # Replace IUPAC codes with corresponding alleles
  func.hmp[func.hmp == "R"] <- "A/G"
  func.hmp[func.hmp == "Y"] <- "C/T"
  func.hmp[func.hmp == "S"] <- "G/C"
  func.hmp[func.hmp == "W"] <- "A/T"
  func.hmp[func.hmp == "K"] <- "G/T"
  func.hmp[func.hmp == "M"] <- "A/C"
  func.hmp[func.hmp == "N"] <- "-"
  
  generate_genefiles <- function(chrom){
    chr <- chrom$Chromosome
    start <- chrom$start
    stop  <- chrom$stop
    gene_id <- chrom$gene_id
    gene_snp <- subset(func.hmp, func.hmp$chrom %in% chr)
    gene_snp <- gene_snp[gene_snp$pos >= start & gene_snp$pos <= stop, c(4, 12:ncol(gene_snp))]
    gene_snp <- t(gene_snp)
    scol <- gene_snp[1, ]
    gene_snp <- data.frame(gene_snp[-1, ])
    colnames(gene_snp) <- scol
    locus <- gene_snp %>% dplyr::select(sort(colnames(gene_snp)))
    write.csv(locus, file = file.path(ann_path, paste0(gene_id, ".csv")))
    return(gene_id)
  }
  
  lapply(1:nrow(g11), function(i) generate_genefiles(g11[i, ]))
  
  return("Processing completed. Check the output directory for results.")
}
