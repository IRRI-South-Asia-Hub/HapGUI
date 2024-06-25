list.of.packages <- (c("profvis","shiny","shinythemes",
                       "shinyFiles","fs", "htmltools","mrMLM",
                       "dplyr","CMplot","tidyr","plyr",
                       "tidyverse","agricolae", "ggplot2",
                       "ggpubr","ggridges","plotly","readxl",
                       "RColorBrewer","crayon","stringr",
                       "berryFunctions","imager","yhat",
                       "corrplot","lme4", "data.table",
                       "ggfortify","extrafont","stats",
                       "psych","plyr","foreach","doParallel","devtools"))

new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

if (!require("petersonR/bestNormalize")){
  devtools::install_github("petersonR/bestNormalize",dependencies = T)
}

if (!require("flextable")){
  install.packages("flextable", type = "binary")
}

list.of.packages <- (c("fastmatch", "haplotypes","augmentedRCBD"))
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages,dependencies = T, 
                                          type = "source")
