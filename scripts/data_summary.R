data_summary <- function(data, varname, groupnames){
  require(plyr)
  
 # groupnames <- as.character(groupnames)
  
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }

  data_sum <- ddply(data, .(get(groupnames)), .fun = summary_func, varname)
  return(data_sum)
}
