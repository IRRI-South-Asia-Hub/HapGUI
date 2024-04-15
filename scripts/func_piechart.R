func_piechart <- function(gene_infile,pheno_file, dir) {

  genes1<- read.csv(gene_infile, header = T)
  colnames(genes1) <- "V1"
  genes <- genes1$V1
  
  all_pheno <- read.csv(pheno_file)

  for (gun in genes) {
    
    colm <- colnames(all_pheno)
    colm<- colm[-1]
    
    for (jun in colm) {

      print(colm)
      haplo_P <- file.path(dir,gun,jun)
      z <- haplo_P
      setwd(paste0(z))
      
      hap_file <- "Haplotype.csv"
      if (file.exists(hap_file)) {
        pie_data <- read.csv(hap_file)
        names(pie_data)[1] <- "Haplotype"
        
        #pie_chart 
        no_haps <- length(unique(pie_data$Haplotype))
        colr <- colorRampPalette(c("orange", "beige", "brown", "gold", "blue", "green"))(no_haps)
        pie_chart1 <- ggplot(pie_data, aes(x = " ", y = Freq, fill = Haplotype,color)) +
          geom_bar(stat = "identity", width = 1) +  scale_fill_manual(values = colr) +
          coord_polar(theta = "y") +
          #scale_fill_manual(values = colr(16)) +
          labs(title = "Haplotype Frequencies", fill = "Haplotypes", x = NULL, y = NULL) +
          theme(plot.title = element_blank(),
                legend.position = "right",
                legend.title = element_text(size = 14),
                panel.background  = element_rect(fill = "white"),
                legend.text = element_text(size = 12),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                plot.margin = margin(1, 1, 1, 1)) +
          geom_text(aes(label = ""),
                    position = position_stack(vjust = 0.5),
                    size = 2,
                    hjust = 0.2,
                    fontface = "bold",
                    family = "Arial",
                    color = "white")
        
        output_file1 <- "hap_diversity.png"
        ggsave(output_file1, pie_chart1, width = 8, height = 6, dpi = 600)
      }
      
      hap_file2 <- "haplotype_subset.csv"
      if (file.exists(hap_file2)) {
        pie_data2 <- read.csv(hap_file2)
        names(pie_data2)[1] <- "Haplotype"
        
        no_haps1 <- length(unique(pie_data2$Haplotype))
        colr1 <- colorRampPalette(c("orange", "beige", "brown", "gold", "blue", "green"))(no_haps1)  
        pie_chart2 <- ggplot(pie_data2, aes(x = "", y = freq_pop, fill = Haplotype)) +
          geom_bar(stat = "identity", width = 1) + scale_fill_manual(values = colr1) +
          coord_polar(theta = "y") +
          #scale_fill_manual(values = colr(9)) +
          labs(title = "Haplotype Frequencies", fill = "Haplotypes", x = NULL, y = NULL) +
          theme(plot.title = element_blank(),
                legend.position = "right",
                legend.title = element_text(size = 14),
                panel.background  = element_rect(fill = "white"),
                legend.text = element_text(size = 12),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                plot.margin = margin(1, 1, 1, 1)) +
          geom_text(aes(label = ""),
                    position = position_stack(vjust = 0.5),
                    size = 2,
                    hjust = 0.2,
                    fontface = "bold",
                    family = "Arial",
                    color = "white")
        
        output_file2 <- "hap_subset_diversity.png"
        ggsave(output_file2, pie_chart2, width = 8, height = 6, dpi = 600)
      }
    }
  }
}
