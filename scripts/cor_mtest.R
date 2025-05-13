cor_mtest <- function(mat, ...) {
  # mat <- as.matrix(mat)
  # n <- ncol(mat)
  # p.mat<- matrix(NA, n, n)
  # diag(p.mat) <- 0
  # for (i in 1:(n - 1)) {
  #   for (j in (i + 1):n) {
  #     tmp <- cor.test(mat[, i], mat[, j], ...)
  #     p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
  #   }
  # }
  # colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  # p.mat

  
  pairs.panels(mat,
               smooth = TRUE,      # If TRUE, draws loess smooths
               scale = FALSE,      # If TRUE, scales the correlation text font
               density = TRUE,     # If TRUE, adds density plots and histograms
               ellipses = TRUE,    # If TRUE, draws ellipses
               method = "pearson", # Correlation method (also "spearman" or "kendall")
               pch = 21,           # pch symbol
               lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
               cor = TRUE,         # If TRUE, reports correlations
               jiggle = FALSE,     # If TRUE, data points are jittered
               factor = 2,         # Jittering factor
               hist.col = 3,       # Histograms color
               stars = TRUE,       # If TRUE, adds significance level with stars
               ci = TRUE)          # If TRUE, adds confidence intervals
  
}