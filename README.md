# HaploGUI
install following libraries 





```R

install.packages(c("mrMLM", "dplyr", "CMplot", "bestNormalize", "tidyr", "plyr", "tidyverse", 
                   "haplotypes", "agricolae", "ggplot2", "ggpubr", "plotly", "readxl", 
                   "RColorBrewer", "crayon", "stringr", "berryFunctions", "imager", 
                   "yhat", "shiny", "shinythemes", "corrplot", "augmentedRCBD", "lme4", 
                   "data.table", "ggfortify", "extrafont", "stats", "psych","BiocManager","rJava","devtools"))
```
After running the following codes we need to run the following command to view the web app
run the following code

```R
shiny::runGitHub(repo = "HaploGUI","IRRI-South-Asia-Hub")
```

If you want to see all the files at your desired directory,u can use 

```R
shiny::runGitHub(repo = "HaploGUI","IRRI-South-Asia-Hub",destdir = "D:/foldername/foldername2")
```
