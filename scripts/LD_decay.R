ld_decay <- function(vcf_path, ip_dir, dir) {
  message("Running LD-decay analysis...")
  
  # Run PopLDdecay command
  # if (choice_geno == "option1") {
  
  system(command = paste0(ip_dir, "/PopLDdecay/bin/PopLDdecay ",
                          "-InVCF ", vcf_path,
                          " -OutStat ", file.path(dir, "marker_LD_stats"),
                          " -MaxDist 500 -OutType 1"))

  ld_file <- file.path(dir, "marker_LD_stats.stat.gz")
  ld_data <- read.table(ld_file, header = FALSE, sep = "\t")
  colnames(ld_data) <- c("Dist", "r2", "Mean_D", "Sum_r^2", "Sum_D", "NumberPairs")
  
  write.csv(ld_data, file.path(dir, "marker_LD_stats.csv"), row.names = FALSE)
  
  p <- ggplot(ld_data, aes(x = Dist, y = r2)) +
    geom_point(color = "white", alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", span = 0.2, color = "black", se = FALSE) +
    labs(
      title = "LD Decay",
      x = "Distance (Kb)",
      y = expression(mean~r^2)
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold", color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      axis.line = element_line(color = "black", linewidth = 0.8)
    )
  
  ggsave(file.path(dir, "LD_plot.png"), plot = p, width = 7, height = 6, dpi = 600)
  
  return(ld_data)
#   }else{
#     
#     system(command = paste0(ip_dir, "/PopLDdecay/bin/PopLDdecay ",
#                             "-InVCF ", vcf_path,
#                             " -OutStat ", file.path(dir, "marker_LD_stats"),
#                             " -MaxDist 500 -OutType 1"))
#     
#     ld_file <- file.path(dir, "marker_LD_stats.stat.gz")
#     ld_data <- read.table(ld_file, header = FALSE, sep = "\t")
#     colnames(ld_data) <- c("Dist", "r2", "Mean_D", "Sum_r^2", "Sum_D", "NumberPairs")
#     
#     write.csv(ld_data, file.path(dir, "marker_LD_stats.csv"), row.names = FALSE)
#     
#     p <- ggplot(ld_data, aes(x = Dist, y = r2)) +
#       geom_point(color = "white", alpha = 0.3, size = 1) +
#       geom_smooth(method = "loess", span = 0.2, color = "black", se = FALSE) +
#       labs(
#         title = "LD Decay",
#         x = "Distance (Kb)",
#         y = expression(mean~r^2)
#       ) +
#       theme_bw(base_size = 14) +
#       theme(
#         plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
#         axis.title = element_text(size = 16, face = "bold", color = "black"),
#         axis.text = element_text(size = 14, color = "black"),
#         panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
#         axis.line = element_line(color = "black", linewidth = 0.8)
#       )
#     
#     ggsave(file.path(dir, "LD_plot.png"), plot = p, width = 7, height = 6, dpi = 600)
#   }
#   
#   
#   
}

# ld_data <- ld_decay(
#   vcf_path = file.path(dir, "marker.vcf"),  
#   ip_dir = ip_dir,                         
#   dir = dir                                 
# )
