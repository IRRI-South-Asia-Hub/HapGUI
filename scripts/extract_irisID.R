extract_irisID <- function(trait, infile, perc, dir, ip_dir){
  
  theme_cus <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                     strip.background=element_blank(),axis.text.x=element_text(colour="black"),
                     axis.text.y=element_text(colour="black"),
                     axis.ticks=element_line(colour="black"),
                     plot.margin=unit(c(1,1,1,1),"line"),
                     legend.title = element_text(color = "black", size = 10),
                     legend.text = element_text(color = "red"),
                     legend.background = element_rect(fill = "lightgray"),
                     legend.key = element_rect(fill = "white", color = NA))
  phenofile <- infile
  # Step1: Phenotypic extract -----------------------------------------------
  
  dir.create(file.path(dir,"EtGWAS_results"))
  data_dir <- paste0(dir,"/EtGWAS_results/")
  
  outlist <- paste0(data_dir,trait,"_list.txt")
  hist_plot <- paste0(data_dir,trait,"_histogram.png")
  box_plot <- paste0(data_dir,trait,"_boxplot.png")
  
  phe <- read.csv(infile, header = T)
  phe <- phe[,1:2]
  colnames(phe) <- c("Accessions","trait")
  # phe$Accessions <- gsub(pattern = "IRIS ",replacement = "IRIS_",
  #                       phe$Accessions)
  
  box <- ggplot(phe) + aes(y = trait) +
    geom_boxplot(fill = "#0c4c8a") + theme_cus + xlab("") + ylab(trait)
  ggsave(box, filename = box_plot,width = 10,height = 6)
  
  his <- ggplot(data = phe, aes(x = trait)) + geom_histogram()+
    theme_cus + xlab("")
  
  ggsave(his,filename = hist_plot,width = 10,height = 6)
  
  out1 <- phe[,c(1,1)]
  write.table(out1, file = outlist, col.names = F, row.names = F,
              quote = F)
  
  
  # # Step 2: genotypic data extract------------------------------------------------
  infam <- "marker_etgwas"
  files_to_copy <- list.files(pattern = "^marker_etgwas", full.names = TRUE)
  file.copy(from = files_to_copy, to = data_dir, overwrite = TRUE)
  
  system(paste0(ip_dir,"/plink2 --bfile ",
                data_dir,infam," --export ped --out ",data_dir,infam))
  
  remove(out)
  
  # pca_plink ---------------------------------------------------------------
  
  theme_cus <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                    strip.background=element_blank(),
                    axis.text.y=element_text(colour="black", face = "bold"),
                    axis.text.x=element_text(colour="black", face = "bold",angle = 90),
                    axis.ticks=element_line(colour="black"),
                    strip.text.x = element_text(colour="black", face = "bold"),
                    title = element_text(color = "black", size = 12,
                                         face = "bold"),
                    plot.margin=unit(c(1,1,1,1),"line"),
                    legend.title = element_text(color = "black", size = 12, 
                                                face = "bold"),
                    legend.text = element_text(colour = "black", size = 12,),
                    legend.key = element_rect(fill = "white", color = NA))
  
  # pheno_dist_xp -----------------------------------------------------------
  # cat(green("Pheno_dist\n"))
  # 
  # popfile <- "Edited all population.txt"
  # 
  # Type <- c("indx","ind2","aus","ind1B","ind1A",
  #           "ind3","japx","temp","trop","subtrop","admix","aro")
  # 
  # Initial distribution ----------------------------------------------------
  b <- read.delim(phenofile, header = T, sep = ",")
  b <- b[,1:2]
  names(b) <- c("Designation","phenotype")
  val <- ceiling((perc*nrow(b))/100)
  
  low <- max(head(b[order(b$phenotype, decreasing=F), ], val)[,2])
  high <- min(head(b[order(b$phenotype, decreasing=T), ], val)[,2])
  
  b$group <- ifelse(b$phenotype <= low, 1, ifelse(b$phenotype >= high, 3, 2))
  supp.labs <- c("Low bulk", "Mixed","High bulk")
  names(supp.labs) <- c("1", "2", "3")
  
  image_file <- paste0(data_dir,trait,"_",perc,"_dist1.png")
  #png(image_file, width = 1000, height = 650)
  dist <- ggplot(data=b, aes(phenotype)) + theme_bw() + 
    theme(legend.position="none",
          panel.background = element_blank(),panel.border=element_rect(fill=NA),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          strip.background=element_blank(),
          panel.spacing = unit(0.3, "lines"),
          strip.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.ticks = element_line(colour = "black"),
          axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(colour="black", size = 12)) +
    geom_histogram(data=subset(b,group==1), fill = "red", color = "red", alpha = 0.5)+
    geom_histogram(data=subset(b,group==2), fill = "green", color = "green", alpha = 0.5) +
    geom_histogram(data=subset(b,group==3), fill = "blue", color = "blue", alpha = 0.5)+
    facet_grid(. ~ group, scales = "free", labeller = labeller(group = supp.labs))
  dist
  #dev.off()
  ggsave(image_file, dist, width = 8, height = 6)
  
  pheno_out <- paste0(data_dir,trait,"_",perc,"_pheno.txt")
  write.table(b, file = pheno_out, sep = "\t", row.names = F, 
              col.names = T, quote = F)
  
  #Extra
  low <- b[b$group==1,c(1,1)]
  write.table(low, file = paste0(data_dir,trait,"_low",perc,".txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  high <- b[b$group==3,c(1,1)]
  write.table(high, file = paste0(data_dir,trait,"_high",perc,".txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  rand <- b[b$group==2,c(1,1)]
  rand2 <- rand[as.integer(runif(min(nrow(low),nrow(high)),1,nrow(rand))),]
  write.table(rand2, file = paste0(data_dir,trait,"_rand",perc,".txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  cat(green("Completed the pheno analysis\n"))
  
  
  # Pooling genotype --------------------------------------------------------
  
  low_in <- paste0(data_dir,trait,"_low",perc,".txt")
  rand_in <- paste0(data_dir,trait,"_rand",perc,".txt")
  high_in <- paste0(data_dir,trait,"_high",perc,".txt")
  
  low_out <- paste0(data_dir,trait,"_low",perc)
  rand_out <- paste0(data_dir,trait,"_rand",perc)
  high_out <- paste0(data_dir,trait,"_high",perc)
  
  system(paste0(ip_dir,"/plink2 --bfile ",data_dir,
                infam," --keep ",low_in," --export ped --out ",low_out))
  
  system(paste0(ip_dir,"/plink2 --bfile ",data_dir,
                infam," --keep ",rand_in," --export ped --out ",rand_out))
  
  system(paste0(ip_dir,"/plink2 --bfile ",data_dir,
                infam," --keep ",high_in," --export ped --out ",high_out))
  
  allele_file <- paste0(data_dir,trait,"_hmp_allele.txt")
  
  system(command = paste0(ip_dir, "/plink2 --bfile ",data_dir, infam," --freq --out et_marker"))
  
  system(paste0("bash scripts/pooling_snp.sh -h ",allele_file))
  
  bulks <- c("high","low","rand")
  
  pooling_snp <- function(i){
    
    infile <- paste0(data_dir,trait,"_",bulks[i],perc,".ped")
    c <- read.delim(infile, header = F, sep = "\t", colClasses = c("character"))
    d <- c[,c(7:ncol(c))]
    out <- apply(d, 2, function(x){
      run = rle(sort(x[x!=0]))
      if(length(run$lengths) < 2){
        run$lengths[2] <- 0
        run$values[2] <- 0
      } 
      g <- stack(run)[1]
    })
    
    lstData <- Map(as.data.frame, out)
    dfrData <- rbindlist(lstData)
    
    count_cols <- sort(c(seq(1,nrow(dfrData),4),seq(2,nrow(dfrData),4)))
    all_cols <- sort(c(seq(3,nrow(dfrData),4),seq(4,nrow(dfrData),4)))
    df <- dfrData[count_cols]
    df$all <- dfrData[all_cols]
    
    colhead <- data.frame(rep(seq(1,(nrow(df)/4)), each=4))
    colnames(colhead) <- "num"
    
    comb <- cbind(colhead,df)
    comb$values <- as.numeric(comb$values)
    comb$alle <- paste(comb$num, comb$all, sep = "_")
    comb <- comb[comb$all != 0,]
    out <- comb %>% dplyr::group_by(alle) %>% 
      dplyr::summarize(tot = sum(values), nums = unique(num), .groups = 'drop')
    
    outfile <- paste0(data_dir,"temp_",bulks[i],".csv")
    write.csv(out, outfile, row.names = F)
  }
  cl<- detectCores()
  registerDoParallel(cl)
  result <- foreach (i=1:3, .packages = c("dplyr","data.table")) %dopar% {
    pooling_snp(i)
  }
  stopImplicitCluster()
  
  print("Finshed the pooling function")
  
  out_low <- read.csv(paste0(data_dir,"temp_low.csv"))
  out_high <- read.csv(paste0(data_dir,"temp_high.csv"))
  out_rand <- read.csv(paste0(data_dir,"temp_rand.csv"))
  
  temp <- merge(out_high,out_low, by = "alle", all = T)
  comb <- merge(temp,out_rand, by = "alle", all = T)
  out <- comb[order(comb$nums.x),c(1,2,4,6,7)]
  colnames(out) <- c("RefAlt", "high" ,"low","random","num")
  
  # Ref and Alt alleles From HapMap -----------------------------------------
  
  hmp <- read.delim(allele_file,header = T,sep = "\t")
  hmp$num <- c(1:nrow(hmp))
  hmp$Ref1 <- paste(hmp$num,hmp$Ref,sep = "_")
  hmp$Alt1 <- paste(hmp$num,hmp$Alt,sep = "_")
  ref <- hmp[,c(3,4)]
  alt <- hmp[,c(3,5)]
  
  ref_comb <- merge(ref,out,by.x = "Ref1", by.y = "RefAlt")
  alt_comb <- merge(alt,out,by.x = "Alt1", by.y = "RefAlt")
  names(alt_comb)[3:5] <- c("high_alt","low_alt","random_alt")
  names(ref_comb)[3:5] <- c("high_ref","low_ref","random_ref")
  
  comb <- merge(ref_comb,alt_comb, by = "num.x", all = T)
  final <- comb[,c(2,7,3,8,4,9,5,10)]
  colnames(final) <- c("Ref","Alt","high_ref","high_alt","low_ref",
                       "low_alt","random_ref","random_alt")
  
  final_out <- paste0(data_dir,trait,"_input",perc,".txt")
  write.table(final, file = final_out, sep = "\t",col.names = T, row.names = F,
              quote = F)
  print("the input file has been created")
  
  # Association -------------------------------------------------------------
  mapfile <- paste0(data_dir,infam,".map")
  
  outfile <- paste(data_dir,trait,"_intermediate_result",perc,".csv",sep = "")
  outqtls <- paste(data_dir,trait,"_qtls",perc,".csv",sep = "")
  outsnps <- paste(data_dir,trait,"_Final_result",perc,".csv",sep = "")
  outpos <- paste(data_dir,"/Et-GWAS","_pos_",trait,perc,".csv",sep = "")
  
  input <- final
  input[is.na(input)]=0
  
  cols = c(3:8)    
  input[,cols] = apply(input[,cols], 2, function(x) as.numeric(as.character(x)));
  map <- read.delim(mapfile, header = F)
  input2 <- cbind(map[,c(2,1,4)],input[,c(3:8)])
  input <- input2
  names(input)[1] <- "snpid"
  names(input)[2] <- "chr"
  names(input)[3] <- "pos"
  
  remove(final)
  xpgwas_modified <- function(input, filter=50,  plotlambda=TRUE){
    statout <- get_chistat(snps=input, filter=filter)
    
    lam <- estlambda(statout$stat, plot = plotlambda, proportion = 1, 
                     method = "regression",
                     filter = TRUE, main="Before genomic control")
    
    outqval <- xpgwas_qval(stat=statout, lambda=lam[['estimate']])
    return(outqval)  
  }
  
  qval <- xpgwas_modified(input, filter=50,  plotlambda=TRUE)
  
  write.csv(qval, file = outfile, row.names = F)
  
  b <- qval[,c(1:3,5)]
  names(b) <- c("SNP","Chromosome","Position","trait1")
  
  suggestiveline = (0.05/nrow(b))*nrow(phe)
  genomewideline = (0.01/nrow(b))*nrow(phe)
  
  CMplot(b,plot.type="m",threshold=c(suggestiveline,genomewideline),
         threshold.col=c('red','orange'),
         multracks=FALSE, chr.den.col=NULL,
         file.name = paste0(trait,"_",perc),
         file="jpg",dpi=600,file.output=TRUE,verbose=TRUE,width=14,height=10)
  
  system(command = paste0("cp *",trait,"_",perc,".jpg ",data_dir,
                          trait,"_",perc,"Manhattan.jpg"))
  
  #xpplot_mine
  
  # Findign significant snps ------------------------------------------------
  
  qval$log10p <- -log10(qval$pval)
  SNPset <- qval
  suggestiveline = -log10((0.05/nrow(b))*nrow(phe))
  
  qtltable <- SNPset[SNPset$log10p >= suggestiveline,]
  
  write.csv(qtltable, file = outsnps, row.names = F)
  
  pos_file <- qtltable[,c(2,3)]
  colnames(pos_file) <- c("Chromosome","Marker.position..bp.")
  
  write.csv(pos_file, file = outpos, row.names = F)
  
  qtltable <- SNPset %>% dplyr::mutate(passThresh = log10p >= suggestiveline) %>%
    dplyr::group_by(chr, run = {
      run = rle(passThresh)
      rep(seq_along(run$lengths), run$lengths)
    }) %>% dplyr::filter(passThresh == TRUE) %>% dplyr::ungroup()
  
  if(nrow(qtltable) > 0){
    qtltable <- qtltable %>% dplyr::ungroup() %>% dplyr::group_by(chr) %>%
      dplyr::group_by(chr, qtl = {
        qtl = rep(seq_along(rle(run)$lengths), rle(run)$lengths)
      }) %>% dplyr::select(-run) %>%
      dplyr::summarize(start = min(pos), end = max(pos), length = end - start + 1,
                       nSNPs = length(pos),
                       avgSNPs_Mb = round(length(pos)/(max(pos) - min(pos)) * 1e+06),
                       .groups = 'drop')
  }
  names(qtltable)[1] <- "CHR"
  write.csv(qtltable, file = outqtls,row.names = F)
  print("everything is done")
  
  list(plot1 = box, plot2 = his, plot5 = dist)
}