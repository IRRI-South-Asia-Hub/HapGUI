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

## 3K PLOTS: The simplest data visualization and analysis tool for 3k population data

### 1. Availability 
The web application named NAME is based on the R shiny package.
It is available on both the shinyapp.io server and the GitHub repository and does not require pre-installation software or packages, except for the GitHub option (see 1.2)

#### 1.1 Shinyapp.io
The web app can be accessed by simply going to the link  (www.link1.com) 

#### 1.2 GitHub
The 3k plots have also been hosted in the GitHub repository of the IRRI South Asia GitHub account.
To be able to run the Shiny app using the hosted GitHub repository we must have the R programming language installed in our device along with a suitable IDE (like Rstudio), further, we need to install the package Shiny using the command install.packages(“shiny”),
More information on how to run the web app is mentioned in section 2.

#### 1.2.1 Installing R
We can go to the Cran.r repository and download the R programming language, later we can go to the posit repository and download the R studio, the following steps can be followed for the same.
https://teacherscollege.screenstepslive.com/a/1108074-install-r-and-r-studio-for-windows

#### 1.2.2 Installing Rshiny

After installation of the R program and R studio, we should open R studio and put the command install.packages(“shiny”)  and run it, This will install the R shiny package.
 




### 2. Pre-requisites 

#### 2.1 Shinyapps.io
To be able to use the web application hosted in the shiny app.io, we only need an internet connection and a browser available on our device.
Simply, click on the link mentioned in section 1 () and we will have the interface right in front of us, read section 3 to understand the components of the interface and section 4 for implementation and data preparation.

2.2 GitHub
To be able to use the web application hosted in the GitHub repository we need an Internet connection, R programming installed, and an IDE(https://teacherscollege.screenstepslive.com/a/1108074-install-r-and-r-studio-for-windows) for R programming.
After all the installation mentioned in section 1.2, we should open the RStudio and use the command runGitHub("3K Plots”, “IRRISAH”).
In the presence of a good internet connection, The above command must open the interface.
3. The Interface
#### 3.1 Sidebar Panels
There are a total of three independent panels with three different sets of functionalities, viz. plots, GWAS and Haplo pheno. We can click on each tab heading to proceed with their functionalities.

#### 3.1.1 Input file for Plots tab
The Input phenotypic file must be in the CSV format (comma delimited values, i.e.,.csv extension), we don’t recommend renaming the .xlsx files into .csv, rather we encourage saving the sheets in a CSV format.

It is necessary to maintain the order of the data in the following order:
1. Location:
The location column should contain the name of the location corresponding to an entry

2. Season:
The location column should contain the name of the Season corresponding to an entry

3. Block:
The name or number representing the block

4. Genotype information column:
Must contain  the genotype information in some form of ID (see below)

5. Traits:

These are the columns where the column header is the name of the trait and column values are numerical values.

The Phenotypic files must contain a Genotypes identification column, which will have accessions for our genotypes, currently, we support three types of IDs:” IRIS_ID”,” IRGC_NO” or ”accession ID “An example file is available at (www.link2.com) 
Note: These IDs must belong to the 3K rice genome (if not it will not generate the topmost plots but will generate the other plots depending on what information you have)
Examples of three IDs:
#### 1.IRIS_ID: IRIS-313-10000
#### 2. IRGC_NO: IRGC 135900
#### 3. ACCESSION: QIUGUANGTENGXI


BLOCKS column
The presence of blocks column and genotypic id plays an important role in the generation of the ANOVA tables and the BLUP values.

Hence we should have a column, that will contain the block information, like block1, block2, or simply 1, 2. etc.

Note: Make sure that the block names are consistent i.e. block1 should be written as block1 everywhere, it should not be like block1 or Block1 or BLOCK1 or block1.

Traits column
This is one of the most important columns, we may have more than one column representing values of traits in proper numeric formats. For generating plots and doing data visualization and analysis we must choose the name of the columns that represent the data for a particular trait. Please refer to the example file at (www.link2.com).

3.1.2 Top Panels
The tab located next to the home tab is the Phenotype Analysis tab.


 
The topmost left section of the plots tab is for uploading the phenotypic file, which is expected to be curated from the user side and the data must be in a similar format as mentioned above in 3.1.1.
In the upload phenotypic file section, we use the browse button to browse files present locally in our device, these files must align with the format mentioned in section 3.1.1.

Select the Genotype ID column section
Just below the browse section, we have radio buttons for three available options representing the forms in which the genotypic information is accepted.
Our data must have at least one column representing the genotype names in terms of Accession or IRIS_ID or IRGC_NO. 
We should choose whatever is available to us, if we have all the IDs, we may choose the column we are most sure about.

Choose genotype column
We need to select the name of the column of our data file which has the genotypic IDs,
For instance if in our data a column is named as IDs and contains all the IRGC Nos then we must select the ID column from the drop-down.

Choose traits column
We need to select the name of the column of our data file which has the trait values, please note we can choose only one trait column at a time.

For instance, in our data, there are 3 columns containing the three phenotypic values (say named: DTF, PH, PL). Then we need to choose any one of them and press generate plots to get all the plots presented on the right side tabs see section 3.1.1.1 for more information.

3.1.1.1 Main panel of plots tab
The main panel of the plots tabs is at the center of the interface and contains several tabs, four of which are for four different data visualization plots, and the other 3 are for summary, ANOVA table, and BLUPs (Adjusted means) respectively.
3.1.1.1.1 Data Visualization tabs
The data visualization tabs display the graphs and also have a download option below, using which we can save the generated plots in a specific directory as we wish.
Note: All the plots will be displayed only after we upload the phenotypic file and fill in all the requirements in the sidebar panel.
3.1.1.1.1.1 Histogram tab
The histogram tab will generate the histogram for the chosen trait, we also have the option to choose colors and then download the same.

3.1.1.1.1.2 Violin tab
The violin plots have three sections, where we can generate three different types of violin plots by adjusting parameters.
1st Violin plot
Given that we have proper genotype IDs, the first plot will generate violin plots with the x-axis separating values based on the 3k subpopulation.

2nd Violin plot
Here we need to choose the name of the column that has data for season or location, if we choose seasons it will give us violin plots based on the seasons present in our data, similarly, it will show for location.
Please Note we should have at least one of the season or location columns if not please ignore the generated violin plot in this section
3rd Violin plot
This is a combined plot where we will get violin plots for both season and location.
To be able to generate season and location-wise plots, we need to choose the location column in the “choose location column” section and choose the season in the “choose season column”
3.1.1.1.1.3 Bar plot tab
Similar to violin plots these also have three sections, we can generate three different types of bar plots by adjusting parameters.
1st bar plot
Given that we have proper genotype IDs, the first plot will generate bar plots with the x-axis separating values based on the 3k subpopulation.
2nd bar plot
Here we need to choose the name of the column that has data for season or location, if we choose seasons it will give us bar plots based on the seasons present in our data, similarly, it will show for location.
Please Note we should have at least one of the season or location columns if not please ignore the generated bar plot in this section
3rd bar plot
This is a combined plot where we will get bar plots for both season and location. To be able to generate season and location-wise plots, we need to choose the location column in the “choose location column” section and choose the season in the “choose season column”
3.1.1.1.1.4 Correlation plot tab
Correlation plots will be generated from the given data for a particular trait, which will make a location season combination corrplot.
To generate a meaningful correlation plot we must choose the season column and the location column using the drop-down.
3.1.1.1.2 Summary tab
The summary table has one row and four columns giving the minimum value, maximum value, mean, median, skewness, kurtosis, and standard deviation of the chosen trait, it can be downloaded as a CSV file.
3.1.1.1.3 ANOVA tab
The ANOVA tab calculates adjusted means using concepts of ANOVA considering the data is in proper augmented RCBD format, moreover, we must specify the column that contains BLOCKS information (see section 4.3). Data in other experimental designs are expected to show errors or produce unreliable results.
The results will be displayed as console output and one can copy the results easily, we don’t recommend using the adjusted means generated from ANOVA hence we use a better model in the next section(3.1.1.1.4)

3.1.1.1.4 BLUPS tab
The BLUPs tab will produce adjusted mean values for all the traits we choose from the checkbox(see section 4). The BLUPS tab has one dropdown where we can choose the column name of our data that contains the block information(block name or number), following this, there is a group of checkboxes where each checkbox represents a column name from our data, we need to choose the traits column names for which we need BLUP values




3.1.2 GWAS tab
The GWAS tabs will accept phenotypic files in some given format(details:3.1.2.1) and generate GWAS results in a folder named results which will be generated in the current working directory.
3.1.2.1 GWAS phenotypic input
The input file for the GWAS must be in CSV format, please refer to section 3.1.1 to understand the CSV format.
The file must contain only two columns where the first one should be named <Phenotype> and the other can be anything.
<Phenotype>
This column must contain the IRIS IDs of the genotypes.
The format should be like “IRIS-313-10000” (Nir)
Column2
The second column must contain the phenotypic values corresponding to the genotypes in the <Phenotype> column, please note, that the program won't entertain missing values.
Following is the example of a sample file:


 
3.1.2.2 Generating results
After we upload our phenotypic files and press the run_gwas button, it will automatically start doing the GWAS analysis. It will take time depending on the no of genotypes.

Two plots and MTAs will be displayed in the respective tabs.
All the results will be produced in the folder named results, which will be generated in our working directory.

3.1.2.2 PCA and Association tabs
The first tab in the main panel of the GWAS tabs(located at the center) is the PCA,
this tab will display the principal components of the given data.
This tab will be followed by Associations tabs which will display the Q-Q plots followed by the Manhattan Plots

3.1.2.2 MTAs tab
The MTA tab will display the Marker Traits Associations generated by the GWAS.
We can also use the download button to save the MTAs locally on our device.

3.1.3 Haplo pheno tab
This tab provides flexibility to the user to use the gene positions generated by the GWAS analysis of the GWAS tab or one may upload the positions in the CSV format which they might have generated by some other means.



3.1.3.1 Input file
1. Position generated from our applications: If we have run the GWAS already the POS.csv file will be generated in our working directory. And we may simply go through the browse button and upload this POS.csv.
2. Positions from the user: The user must have the positions in a CSV format and with two columns, the first column will contain the chromosome numbers and the second will have Marker positions.
3.1.3.3 Candidate genes
As soon as we upload the position files the candidate genes table will be displayed and we can download the generated table in CSV format.
3.1.3.3 Superior haplotypes
LD regions
We can set the LD region for candidate genes using the input box.
High low
We can choose whether or not we want to run haplotype-phenotype analysis for high phenotypic values or low phenotypic values using the radio button.

Run
We can now press the run_hap button which will start the haplo-pheno analysis.
It will generate the result folders in our working directory.

4. Saving results 
Everything can be downloaded using the download buttons present in each section.







 





