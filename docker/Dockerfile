
# Use the latest Rocker Shiny image
FROM rocker/shiny

LABEL maintainer="Bandana Pahi <pahibandana207@gmail.com>"
LABEL author="Bandana Pahi"
LABEL version="1.0"
LABEL description="Created by Bandana Pahi and pushed to irrisah2012's Docker Hub"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    sudo \
    wget \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    software-properties-common \
    libfftw3-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    openjdk-17-jdk \
    default-jdk \
    default-jre \
    ant \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgmp-dev \
    libmpfr-dev \
    libmpc-dev \
    cmake \
    && apt-get clean

# Set up R environment variables
ENV R_LIBS_USER=/usr/local/lib/R/site-library
# Install R dependencies
RUN R -e "install.packages('ragg', repos='https://cloud.r-project.org/', dependencies=TRUE, INSTALL_opts = c('--no-lock'))"
RUN R -e "library(ragg)"
RUN R -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos='https://cloud.r-project.org/')"
RUN R -e "library(remotes)"
RUN R -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools', repos='https://cloud.r-project.org/')"
RUN R -e "library(devtools)"
RUN R -e "install.packages('tidyverse', repos='https://cloud.r-project.org/', dependencies=TRUE, INSTALL_opts = c('--no-lock'))"
RUN R -e "library(tidyverse)"
# Install R packages
RUN R -e "install.packages(c('profvis', 'shiny', 'shinythemes', 'shinyFiles', 'fs', \
    'htmltools', 'mrMLM', 'plyr', 'dplyr', 'CMplot', 'tidyr', 'haplotypes', \
    'agricolae', 'ggplot2', 'ggridges', 'plotly', 'readxl', 'RColorBrewer', \
    'crayon', 'stringr', 'berryFunctions', 'imager', 'yhat', 'corrplot', 'augmentedRCBD', \
    'data.table', 'ggfortify', 'extrafont', 'stats', 'psych', 'foreach', \
    'doParallel', 'qqman', 'reshape2', 'vcfR', 'rJava', 'htmlwidgets'), \
    repos='https://cloud.r-project.org/', dependencies=TRUE, INSTALL_opts = c('--no-lock'))"
# Validate installation
RUN R -e "sapply(c('profvis', 'shiny', 'shinythemes', 'shinyFiles', 'fs', \
    'htmltools', 'mrMLM', 'plyr', 'dplyr', 'CMplot', 'tidyr', 'haplotypes', \
    'agricolae', 'ggplot2', 'ggridges', 'plotly', 'readxl', 'RColorBrewer', \
    'crayon', 'stringr', 'berryFunctions', 'imager', 'yhat', 'corrplot', 'augmentedRCBD', \
    'data.table', 'ggfortify', 'extrafont', 'stats', 'psych', 'foreach', \
    'doParallel', 'qqman', 'reshape2', 'vcfR', 'rJava', 'htmlwidgets'), function(pkg) { \
        if (!requireNamespace(pkg, quietly = TRUE)) \
            message(paste('Package', pkg, 'not found')); \
        })"
RUN R -e "install.packages('bestNormalize', repos='https://cloud.r-project.org/', dependencies=TRUE, INSTALL_opts = c('--no-lock'))"
RUN R -e "if (!requireNamespace('bestNormalize', quietly = TRUE)) remotes::install_github('petersonR/bestNormalize')"
RUN R -e "library(bestNormalize)"
RUN R -e "remotes::install_github('davidgohel/flextable')"
RUN R -e "library(flextable)"
# Install Bioconductor dependencies
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('Rhtslib', 'rtracklayer', 'Rsamtools', 'GenomicAlignments', \
    'BSgenome', 'GenomicFeatures', 'txdbmaker', 'GenomicRanges', 'GenomeInfoDb', \
    'BiocGenerics', 'S4Vectors', 'AnnotationDbi', 'VariantAnnotation'), force = TRUE)"
RUN R -e "sapply(c('Rhtslib', 'rtracklayer', 'Rsamtools', 'GenomicAlignments', \
    'BSgenome', 'GenomicFeatures', 'txdbmaker', 'GenomicRanges', 'GenomeInfoDb', \
    'BiocGenerics', 'S4Vectors', 'AnnotationDbi', 'VariantAnnotation'), \
    function(pkg) { if (!requireNamespace(pkg, quietly = TRUE)) message(paste('Package', pkg, 'not found')) })"
# Install Java and dependencies
RUN apt-get update && apt-get install -y openjdk-17-jdk liblzma-dev libbz2-dev zlib1g-dev
# Set JAVA_HOME
#RUN echo "JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64" >> /etc/environment && \
#    echo "export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64" >> ~/.bashrc && \
#    export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
# Install rJava
RUN R -e "install.packages('rJava', repos='https://cloud.r-project.org/', dependencies=TRUE)"
RUN R -e "Sys.getenv('JAVA_HOME')"
RUN R -e "library(rJava)"
# Install rTASSEL
RUN R -e "if (!requireNamespace('rTASSEL', quietly = TRUE)) { \
              remotes::install_bitbucket('bucklerlab/rTASSEL', ref='master', \
              build_vignettes=FALSE, INSTALL_opts=c('--no-multiarch', '--no-test-load')) }"
#RUN R -e "if (!requireNamespace('rTASSEL', quietly = TRUE)) { \
#              stop('rTASSEL installation failed') }"
RUN R -e "library(rTASSEL)"
RUN apt-get update && apt-get install -y \
    libgfortran5 \
    liblapack-dev \
    libblas-dev
RUN R -e "install.packages(c('Matrix', 'nlme', 'nloptr'), repos='https://cloud.r-project.org/')"
RUN R -e "library(Matrix)"
RUN R -e "library(nlme)"
RUN R -e "install.packages('lme4',repos='https://cloud.r-project.org/', dependencies=TRUE)"
RUN R -e "library(lme4)"
#RUN R -e "install.packages('ggpubr', repos = c("https://cran.rediris.org/", "https://cloud.r-project.org/"), dependencies=TRUE)"
RUN R -e "devtools::install_github('kassambara/ggpubr')"
RUN R -e "library(ggpubr)"
# Create the working directory
RUN mkdir -p /home/app
# Copy Shiny app files into the container
COPY HapGUI /home/app
# Set the working directory
WORKDIR home/app
# Expose the port your Shiny app will run on
EXPOSE 3232
# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/home/app', host='0.0.0.0', port=3232)"]
