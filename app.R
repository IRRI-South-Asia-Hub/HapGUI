setwd("~/Documents/HaploGUI_part2/")
dir <- getwd()


#all libraries needed
library(mrMLM)
library(dplyr)
library(CMplot)
library(bestNormalize)
library(tidyr)
library(plyr)
library(tidyverse)
library(haplotypes)
library(agricolae)
library(ggplot2)
library(ggpubr)
library(plotly)
library(readxl)
library(RColorBrewer)
library(crayon)
library(stringr)
library(berryFunctions)
library(imager)
library(yhat)
library(shiny)
library(shinythemes)
library(corrplot)
library(augmentedRCBD)
library(lme4)
library(data.table)
library(ggfortify)
library(extrafont)
library(stats)
library(psych)
library(plyr)


genes_file <- "ricegenes.txt"

source("data_summary.R")
source("candidate_gene.R")
source("hap_phe.R")
source("func_piechart.R")
source("pca_plot.R")

theme_cus <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                  strip.background=element_blank(),
                  axis.text.x= element_text(size = 14,hjust = 1,angle = 60, colour="black",family = "Arial"),
                  axis.text.y = element_text(size = 14,hjust = 1,colour="black",family = "Arial"),
                  axis.title = element_text(size = 12,face = "bold",family = "Roboto"),
                  legend.text = element_text(size = 12,face="bold",family = "Arial"),
                  legend.title = element_text(color = "black", size = 12),
                  legend.key = element_rect(fill = "white", color = NA),
                  axis.ticks=element_line(colour="black"),
                  plot.margin=unit(c(1,1,1,1),"line"))

mysummary <- function(x,na.rm=F){
  res <- list(Minimum=min(x, na.rm=na.rm),
              Maximum=max(x, na.rm=na.rm),
              Mean=mean(x, na.rm=na.rm),
              StdDev=sd(x,na.rm=na.rm),
              Sknewness = skewness(x),
              Kurtosis = kurtosis(x))
  
  res$Range <- paste0(res$Minimum,"-",res$Maximum)
  res
}

ui = fluidPage(tagList(
  #shinythemes::themeSelector(),
  navbarPage(
    theme = shinythemes::shinytheme("cosmo"),
    id = "my-page",
    title = actionLink("title", "3K-RG Association"),
    tabPanel(
      "Home",icon = icon("home"),
      titlePanel( div(style = "text-align: center;",
                      img(src = "3K PlOTS.png", height = "180px", width = "1200px")))
    ),
    
    #---plots and BLUPs----------
    tabPanel("Phenotype Analysis",fluid = T,
             
             sidebarLayout(
               sidebarPanel(
                 
                 fileInput("pheno","Upload phenotypic data",accept = ".csv"),
                 #Choosing genotype ID type
                 radioButtons("type", "Select genotype_ID Type",selected = "",
                              choices = c("Accession", "IRIS_ID","IRGC_NO")),
                 
                 selectInput("Acc","Choose genotype column",choices = NULL),
                 actionButton("plt_bn","Generate plots"),
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("histogram plot",
                            selectInput("trait","Choose trait column",choices = NULL),
                            tags$div(selectInput("hs_location", "Location", choices = NULL),
                                     style="display:inline-block"),
                            tags$div(selectInput("hs_season", "Season", choices = NULL),
                                     style="display:inline-block"),
                            tags$div(selectInput("hist_cl","Histogram Color",choices =  c(
                              "LightRed" = "#ea5545", "DarkPink" = "#f46a9b",
                              "Orange" = "#ef9b20", "Yellow" = "#edbf33",
                              "LightYellow" = "#ede15b", "Olive" = "#bdcf32",
                              "Green" = "#87bc45", "Blue" = "#27aeef",
                              "Purple" = "#b33dc6", "Red" = "#b30000",
                              "DarkPurple" = "#7c1158", "Violet" = "#4421af",
                              "DeepBlue" = "#1a53ff", "SkyBlue" = "#0d88e6",
                              "Turquoise" = "#00b7c7", "LightGreen" = "#5ad45a",
                              "BrightGreen" = "#8be04e", "Beige" = "#ebdc78",
                              "Salmon" = "#fd7f6f", "greishBlue" = "#7eb0d5",
                              "Lime" = "#b2e061", "DarkPink" = "#bd7ebe",
                              "Apricot" = "#ffb55a", "LightYellow2" = "#ffee65",
                              "Lavender" = "#beb9db", "Pink" = "#fdcce5",
                              "Teal" = "#8bd3c7"
                            ),selected = "#27aeef"),
                            style="display:inline-block"),
                            
                            plotOutput("hist_pl"),
                            downloadButton("hist_dw","Download_histogram")
                   ),
                   tabPanel("violin plot",
                            plotOutput("violin_pl"),
                            downloadButton("violin_dw","Download_violin"),
                            
                            #pooled
                            tags$div(selectInput("location", "Location", choices = NULL),
                                     style="display:inline-block"),
                            tags$div(selectInput("season", "Season", choices = NULL),
                                     style="display:inline-block"),
                            plotOutput("other_violin"),
                            downloadButton("other_violin_dw","Download_violin"),
                   ),
                   
                   tabPanel("bar plot",
                            plotOutput("bar_pl"),
                            downloadButton("bar_dw","Download_bar"),
                            
                            tags$div(selectInput("bar_location", "Location", choices = NULL),
                                     style="display:inline-block"),
                            tags$div(selectInput("bar_season", "Season", choices = NULL),
                                     style="display:inline-block"),
                            plotOutput("other_bar"),
                            downloadButton("other_bar_dw","Download_bar"),
                            
                   ),
                   
                   tabPanel("Corrplot",
                            tags$div(selectInput("corr_location", "Select a location", choices = NULL),
                                     style="display:inline-block"),
                            tags$div(selectInput("corr_season", "Select a season", choices = NULL),
                                     style="display:inline-block"),
                            p("Please select atleast 2 traits"),
                            checkboxGroupInput("corr_traits","Select traits", choices = NULL),
                            plotOutput("corr"),
                            downloadButton("cor_dw","dowload correlation plot")
                   ),
                   
                   tabPanel("Descriptive Statsitics",
                            fluidRow(h3("Overall summary")),
                            tableOutput("table"),
                            downloadButton("summary_dw","download"),
                            fluidRow(h3("Location/Season wise summary")),
                            tableOutput("breakup_table"),
                            downloadButton("break_summary_dw","download")
                   ),
                   
                   tabPanel("ANOVA",
                            h3("ANOVA TABLE"),
                            verbatimTextOutput("anova"),
                            h4("Conclusion"),
                            verbatimTextOutput("con"),
                            downloadButton("aov"),
                            h3("ANOVA for augmented RCBD"),
                            fileInput("nd","augmneted data",accept = ".csv"),
                            actionButton("aug","Generate ANOVA for aug"),
                            selectInput("block","Choose blocks column",choices = NULL),
                            selectInput("id","Treatments",choices = NULL),
                            selectInput("trait2","select trait",choices = NULL),
                            verbatimTextOutput("aug_anova"),
                            tableOutput("means")
                   ),
                   
                   tabPanel("ANOVA new",
                            h3("new anova"),
                            selectInput("anova_treatment","select treatment column",choices = NULL),
                            selectInput("anova_block","select block column",choices = NULL),
                            selectInput("anova_trait","select trait column",choices = NULL),
                            actionButton("runova","run_anova"),
                            verbatimTextOutput("anova_new"),
                   ),
                   
                   tabPanel("BLUPs",
                            tableOutput("blup_output"),
                            downloadButton("blup_dw","download data with ajusted means"),
                            selectInput("block2","Choose blocks column",choices = NULL),
                            checkboxGroupInput("non_trait", "Select traits", choices = NULL),
                            verbatimTextOutput("blupv"),
                            actionButton("genb","Generate BLUPs")
                   )
                 )
               )
             )
    ),
    
    
    # GWAS tab ----------------------------------------------------------------
    tabPanel("GWAS",fluid = T,
             sidebarLayout(
               sidebarPanel(
                 fileInput("gwasfile","Phenotypic file input:"),
                 fileInput("genofile","Genotypic file input:"),
                 fileInput("pcafile","PCA file input:"),
                 actionButton("runa","Run GWAS")
               ),
               mainPanel(
                 tabsetPanel(
                   #------PCA---------------
                   tabPanel("Population structure", 
                            plotOutput("pca_plot", width = "100%", 
                                       height = "300px")),
                   
                   #------Association---------------
                   tabPanel("Association",
                            plotOutput("man"),
                            plotOutput("QQ")
                   ),
                   
                   #----------MTAs------
                   tabPanel("MTAs",fluid = T,
                            tableOutput("result_table")
                   ),
                 )
               )
             )
    ),
    
    
    # hap-phe tab -------------------------------------------------------------
    tabPanel("Haplo-Pheno",fluid = T,
             sidebarLayout(
               sidebarPanel(
                 fileInput("hapfile","Phenotypic file input:"),
                 numericInput("LD", "LD region", 25000,min = 1000, max = 250000),
                 radioButtons("hl","Choose high/low value haplotypes",choices = c("High","Low")),
                 actionButton("runhap","Run Haplo-Pheno")
               ),
               mainPanel(
                 tabsetPanel(
                   #------Candidate genes---------------
                   tabPanel("Candidate Genes",
                            tableOutput("candtable"),
                            downloadButton("canddw")
                   ),
                   tabPanel("Superior haplotypes",
                            tableOutput("haptab"),
                            downloadButton("hapdw","download_haplotype_table"),
                            selectInput("loc_id","Choose Locus ID",choices = NULL),
                            selectInput("season_pie","Choose the season",choices = NULL),
                            imageOutput("piechart")
                   ),
                 )
               )
             )
    )
  )
))

server = function(input, output, session) {
  
  pop =  reactive({
    pop <- read.delim("Edited all population.txt", header = T)
    names(pop) <- c("Designation","Name","Subpopulation","COUNTRY","IRGC")
    return(pop)
  })
  
  #Taking phenotypic file as reactive input
  trait = reactive({
    req(input$pheno)
    a <- read.csv(input$pheno$datapath, header = T)
    return(a)
  })
  
  aug_data =  reactive({
    req(input$nd)
    a<-read.csv(input$nd$datapath,header = T)
    return(a)
  })
  
  #updating the choices 
  observeEvent(c(input$type),{
    locations <- trait()[,1]
    seasons <- trait()[,2]
    
    updateSelectInput(session,"Acc",choices = names(trait()))

    #Histogram
    updateSelectInput(session,"hs_location",choices = locations)
    updateSelectInput(session,"hs_season",choices = seasons)
    updateSelectInput(session,"trait",choices = names(trait())[5:ncol(trait())])
    
    
    updateSelectInput(session,"block2",choices = names(trait()))
    updateCheckboxGroupInput(session,"non_trait",choices = names(trait()))
    updateSelectInput(session,"location",choices = names(trait()))
    
    #violin LS
    updateSelectInput(session,"season",choices = append(names(trait()),"null",0))
    updateSelectInput(session,"location",choices = append(names(trait()),"null"))
    
    #bar LS
    updateSelectInput(session,"bar_season",choices = append(names(trait()),"null",0))
    updateSelectInput(session,"bar_location",choices = append(names(trait()),"null"))
    
    #corr
    updateSelectInput(session,"corr_location",choices = locations)
    updateSelectInput(session,"corr_season",choices = seasons)
    updateCheckboxGroupInput(session,"corr_traits",choices = names(trait())[5:ncol(trait())])
    
    
    #new anova update
    updateSelectInput(session,"anova_treatment",choices = names(trait()))
    updateSelectInput(session,"anova_block",choices = names(trait()))
    updateSelectInput(session,"anova_trait",choices = names(trait()))
  })
  
  #updating the choices for ANOVA
  observeEvent(input$nd,{
    updateSelectInput(session,"block",choices = names(aug_data()))
    updateSelectInput(session,"trait2",choices = names(aug_data()))
    updateSelectInput(session,"id",choices = names(aug_data()))
  })
  
  #Histo--------------------
  #updating the choices for genotypic ID and make single pheno
  data = reactive({
    req(input$pheno)
    trait_d = trait()
    names(trait_d)[which(names(trait_d)== input$trait)] <- c ("phenotype")
    names(trait_d)[which(names(trait_d)== input$Acc)] <- c ("Name")
    
    ycol <- "Name"
    
    if(input$type == "IRIS_ID"){
      ycol <- "Designation" 
    }else if(input$type == "IRGC_NO"){
      ycol <- "IRGC" 
    }
    
    aa <- merge(trait_d, pop(), by.x = "Name", by.y = ycol, all = F)
    return(aa)
  })
  
  histo =  reactive({
    pheno_d2 <- data()[(data()[,2] == input$hs_location & 
                          data()[,3] == input$hs_season),]
    
    p1 <- ggplot(data= pheno_d2, aes(phenotype)) + theme_bw() +
      theme(axis.line = element_line(size=1, colour = "black"),
            panel.grid = element_blank(), panel.border = element_blank(), 
            panel.background = element_blank(), panel.spacing = unit(0.3, "lines"),
            axis.text = element_text(colour="black", size = 12,family = "Arial"),
            #axis.text.y=element_text(colour="black", size = 12),
            axis.title=element_text(size=12,face = "bold",color ="black"),
            strip.text.x = element_text(size = 12, face="bold"))+
      geom_histogram(data=pheno_d2,color="black", fill=input$hist_cl, 
                     alpha = 0.5)  + xlab(input$trait) + ylab("Count")
    return(p1)
  })
  
  output$hist_pl = renderPlot({
    req(input$plt_bn,input$hist_cl)
    histo()
  })
  
  output$hist_dw = downloadHandler(
    filename = function(){
      paste0(input$trait,input$hs_location,input$hs_seasons," histrogram.png")
    },
    content = function(file){
      ggsave(file,histo(), width = 8, height = 6,dpi = 600)
    }
  )
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #Pheno prep with subpop----
  ##updating the choices for genotypic ID and make single pheno with subpops 
  data_sub =  reactive({
    
    req(input$pheno)
    as <- read.csv(input$pheno$datapath)
    names(as)[which(names(as)== input$trait)] <- c ("phenotype")
    names(as)[which(names(as)== input$Acc)] <- c ("Name")
    ycol <- "Name"
    
    #This allows people to choose whether they have IRS_ID or Accession file
    if(input$type == "IRIS_ID"){
      ycol <- "Designation" 
    }else if(input$type == "IRGC_NO"){
      ycol <- "IRGC" 
    }
    
    aas <- merge(as, pop(), by.x = "Name", by.y = ycol, all = F)
    
    subpop <- c("indx", "ind2", "ind1B", "ind1A", "ind3", "aus", "japx", "temp", "trop", "subtrop", "admix", "aro")
    subcol <- c("ind", "ind", "ind", "ind", "ind", "aus", "jap", "jap", "jap", "jap", "admix", "aro")
    aas$sub <- aas$Subpopulation
    for (i in c(1:12)) {
      aas$sub[aas$Subpopulation == subpop[i]] <- subcol[i]
    }
    return(aas)
  })
  
  sub_color = reactive({
    custom_colors = c("jap" = "#E6B0AA","aus" = "#D7BDE2" , "ind" = "#A9CCE3","admix" = "#CCD1D1" ,"aro" = "#A9DFBF")
  })
  
  subpopulation_color = reactive({
    subpop_clr = c( "temp" = "#F1948A","aus" ="#D7BDE2","trop"="#E6B0AA",  "indx"= "#D4E6F1","ind1A"= "#2E86C1","ind1B"= "#5499C7",   "admix"= "#CCD1D1",   "ind2"= "#85C1E9",    "subtrop"="#CD6155", "japx"="#E74C3C","ind3"= "#A9CCE3",   
                    "aro" = "#A9DFBF" )
  })
  
  #voilin tab ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  violin = reactive({
    p2 <- ggplot(data_sub(), aes(x = sub,y=phenotype, fill=sub)) + 
      geom_violin()+ theme_cus +
      #scale_fill_discrete(name = "Subpopulation")
      ylab(input$trait) + geom_boxplot(width=0.1) + 
      scale_fill_manual("",values = sub_color() )+xlab(" ")
    return(p2)
  })
  
  violin_others <- reactive({
    p2 <- ggplot(data_sub(), aes(x = data_sub()[[input$location]],y=phenotype, fill= data_sub()[[input$season]])) + 
      geom_violin()+ theme_cus +
      # scale_fill_discrete(name = input$other_cols)
      ylab(input$trait) + geom_boxplot(width=0.1,position = position_dodge(0.9)) + scale_fill_discrete(name = input$season)+xlab(" ")
    #scale_fill_manual("",values = sub_color() )+xlab(" ")
    return(p2)
  })
  
  output$violin_pl = renderPlot({
    req(input$plt_bn)
    violin()
  })
  
  output$other_violin = renderPlot({
    req(input$plt_bn,input$location, input$season)
    violin_others()
  })
  
  output$violin_dw = downloadHandler(
    filename = function(){
      paste0(input$trait," violin_plot.png")
    },
    content = function(file){
      ggsave(file,violin(), width = 8, height = 6,dpi = 600)
    }
  )
  
  output$other_violin_dw = downloadHandler(
    filename = function(){
      paste0(input$trait,"_",input$location,input$season," violin_plot.png")
    },
    content = function(file){
      ggsave(file,violin_others(), width = 8, height = 6,dpi = 600)
    }
  )
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #bar tab----
  bar = reactive({
    df <- data_summary(data_sub(), varname="phenotype", 
                       groupnames=c("Subpopulation"))
    df <- rename(df, c("phenotype" = "mean"))
    
    p3 <- ggplot(df, aes(x = Subpopulation,y=phenotype, fill=Subpopulation)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=phenotype-sd, ymax=phenotype+sd), width=.2,
                    position=position_dodge(.9)) + theme_cus + ylab(input$trait)+
      scale_fill_manual("",values = subpopulation_color() )
    return(p3)
  })
  
  bar_others = reactive({
    df <- data_summary(data_sub(), varname="phenotype", 
                       groupnames=c(input$bar_location,input$bar_season))
    
    df <- rename(df, c("phenotype" = "mean"))
    p3 <- ggplot(df, aes(x = df[[input$bar_location]],y=phenotype, fill=df[[input$bar_season]])) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=phenotype-sd, ymax=phenotype+sd), width=.2,
                    position=position_dodge(.9)) + theme_cus + ylab(input$trait)+
      scale_fill_discrete(name = df[[input$bar_season]]) +xlab("")+
      labs(fill = "Season", color="Season")
    #scale_fill_manual("",values = subpopulation_color() )
    return(p3)
  })
  
  output$bar_pl = renderPlot({
    req(input$plt_bn)
    
    bar()
    
  })
  
  output$other_bar <- renderPlot({
    req(input$plt_bn,input$bar_location, input$bar_season)
    bar_others()
  })
  
  output$bar_dw = downloadHandler(
    filename = function(){
      paste0(input$trait,"_barplot.png")
    },
    content = function(file){
      ggsave(file,bar(), width = 8, height = 6, dpi = 600)
    }
  )
  
  output$other_bar_dw = downloadHandler(
    filename = function(){
      paste0(input$trait,"_",input$bar_location,input$bar_season,"_barplot.png")
    },
    content = function(file){
      ggsave(file,bar_others(), width = 8, height = 6, dpi = 600)
    }
  )
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #correlation----
  corr = reactive({
    aa <- trait()[(trait()[,1] == input$corr_location & 
                     trait()[,2] == input$corr_season),]
    
    dd <- aa[,input$corr_traits]
    print(dd)
    
    png(file.path(dir,paste0(input$corr_location,input$corr_season,"_corr.png")), 
        width = 6, height = 6, units = "in", res = 300)
    
    pairs.panels(dd, smooth = TRUE, scale = FALSE,
                 density = TRUE, ellipses = TRUE, method = "pearson", 
                 pch = 21, lm = FALSE, cor = TRUE, jiggle = FALSE,
                 factor = 2, hist.col = 3, stars = TRUE, ci = TRUE)       
    dev.off()
  })
  
  output$corr <- renderImage({
    outfile <- tempfile(fileext = '.png')
    corr()
    # Return a list containing the filename
    list(src = file.path(dir,paste0(input$corr_location,input$corr_season,
                                    "_corr.png")), 
         contentType = 'image/png',
         width = 400,
         height = 400,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  output$cor_dw = downloadHandler(
    filename <- function() {
      paste0(input$corr_location,input$corr_season,"_corr.png")
    },
    
    content <- function(file) {
      file.copy(file.path(dir,paste0(input$corr_location,input$corr_season,
                                     "_corr.png")), file)
    },
    contentType = "image/png"
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  summ = reactive({
    
    column_data <- trait()[,c(4:ncol(trait()))]
    colnames(column_data)[1] <- "Name"
    summary_stats <- melt(setDT(column_data),
                          id.vars = c("Name"))[,mysummary(value),
                                               by=.(variable)]
    summary_stats <- rename(summary_stats,"Trait"=variable)
    
    return(summary_stats)
  })
  
  break_summ = reactive({
    column_data <- trait()[,-c(3,4)]
    summary_stats <- melt(setDT(column_data),
                          id.vars=c("LOCATION","SEASON"))[,mysummary(value),
                                                          by=.(variable,LOCATION,SEASON)]
    summary_stats <- rename(summary_stats,"Trait"=variable)
    return(summary_stats)
  })
  
  output$table = renderTable({
    summ()
  })
  
  output$breakup_table = renderTable({
    break_summ()
  })
  
  output$summary_dw = downloadHandler(
    filename = function(){
      "summary.csv"
    },
    content = function(file){
      write.csv(summ(),file)
    }
  )
  
  output$break_summary_dw = downloadHandler(
    filename = function(){
      "loc_season_summary.csv"
    },
    content = function(file){
      write.csv(break_summ(),file)
    }
  )
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #aug anova-------------------- 
  anova_tab = reactive({
    ds = data()
    ds = select_if(ds,is.numeric)
    comb  = gather(data = ds,key = "Factors",value = "Values")
    results = anova(aov(Values~Factors,data = comb))
    res.aov3 <- aov(len ~ supp + dose + supp:dose, data = my_data)
    return(results)
  })
  
  output$anova =  renderPrint({
    req(input$pheno,input$Acc,input$trait)
    anova_tab()
  })#include.rownames = TRUE)
  
  #Basic ANOVA conclusions
  output$con = renderText({
    req(input$pheno,input$Acc,input$trait)
    results = anova_tab()
    
    if(results$`Pr(>F)`[1]<0.05){
      print("the differences in factors are significant")
    }else{
      print("the differences are not significant")
    }
  })
  
  #Augmented ANOVA    
  output$aug_anova = renderPrint({
    d = aug_data()
    d = na.omit(d)
    print(class(d))
    ph =  input$trait2
    Block = input$block
    id = input$id
    DAU.test(d$Block,d$id,d$ph,method = "lsd",console = T,group =F)
  })
  
  output$means = renderTable({
    d = read.delim(input$nd$datapath)
    d = na.omit(d)
    ph =  d$input$trait2
    Block = d$input$block
    id = d$input$id
    out = DAU.test(block = Block,trt = id,y = ph,method = "lsd",console = F,group =F)
    return(out$mean)
  })
  
  #NEW anova
  anova_overall = reactive({
    ds = trait()
    ds = na.omit(ds)
    
    Treatment = ds[[input$anova_treatment]]
    Block = ds[[input$anova_block]]
    trait = ds[[input$anova_trait]]
    
    Treatment = as.factor(Treatment)
    Block = as.factor(Block)
    
    augrc = augmentedRCBD(Block,Treatment,trait,method.comp = "lsd",
                          alpha = 0.05, group = FALSE , console = TRUE)
    return(augrc)
  })
  
  output$anova_new = renderPrint({
    req(input$runova)
    return(anova_overall()$'ANOVA')
  })
  
  output$aov = downloadHandler(
    filename = function(){
      "anova.csv"
    },
    content = function(file){
      write.csv(anova_tab(),file)
    }
  )
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #BLUPs----
  output$blup_output = renderTable({
    #introducing req button to reduce error
    req(input$genb)
    
    d = trait()
    #removing nas to keep balance 
    d = na.omit(d)
    
    #making blocks as factors
    d[[input$block2]] = as.factor(d[[input$block2]])
    
    acc_group <- d[[input$Acc]]
    block = d[[input$block2]]
    
    l = list()
    for (i in input$non_trait){
      model<-lmer(d[[i]]~(1|acc_group)+ (1|block),data = d)
      blup_val <- as.data.frame(coef(model)$acc_group)
      names(blup_val)[which(names(blup_val)=="(Intercept)")] = paste0(i,"(adjusted)")
      
      #list making for column bind
      l[[i]] = blup_val
      #print(l)
    }
    
    #creating a good d
    k = bind_cols(l)
    k$"Genotype" = row.names(k)
    final =  merge(d,k,by.x = input$Acc,by.y = "Genotype")
    return(head(final))
  })
  
  #introducing BLUP reactive for easy downloads
  blup_out = reactive({
    #introducing req button to reduce error
    req(input$genb)
    d = trait()
    #removing nas to keep balance 
    d = na.omit(d)
    #making blocks as factors
    d[[input$block2]] = as.factor(d[[input$block2]])
    acc_group <- d[[input$Acc]]
    block = d[[input$block2]]
    
    l = list()
    for (i in input$non_trait){
      model<-lmer(d[[i]]~(1|acc_group)+(1|block),data = d)
      blup_val <- as.data.frame(coef(model)$acc_group)
      names(blup_val)[which(names(blup_val)=="(Intercept)")] = paste0(i,"(adjusted)")
      #list making for column bind
      l[[i]] = blup_val
    }
    
    #creating a good d
    k = bind_cols(l)
    k$"Genotype" = row.names(k)
    
    final =  merge(d,k,by.x = input$Acc,by.y = "Genotype")
    return(final)
  })
  
  output$blup_dw = downloadHandler(
    filename = function(){
      "blups.csv"
    },
    content = function(file){
      write.csv(blup_out(),file)
    }
  )
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #GWAS tabs----
  phe = reactive({
    phe <- read.csv(file = input$gwasfile$datapath,header = T) #file name  =  pheno.csv
    
    pnam <- c("<Phenotype>","trait") #Change trait
    colnames(phe) <- pnam
    phe$"<Phenotype>" <- gsub(pattern = "IRIS ",replacement = "IRIS-",
                           phe$"<Phenotype>")
    return(phe)
  })
  
  output$analysis_started = renderPrint({
    req(input$runa)
    cat("analysis started","\n")
  })
  
  plotData <- reactiveVal(NULL)
  
  observeEvent(input$runa,{
    req(input$genofile)
    req(input$pcafile)
    
    showModal( modalDialog(
      h4(paste0("GWAS for ",nrow(phe())," genotypes:")),
      footer=tagList(h3("running..."))
    ))
    
    dir.create(file.path(dir,"results"))
    write.csv(phe() ,file.path(dir,"results",file = "pheno.csv"),row.names = F)
    
    pca_plot(input$pcafile$datapath) 
    
    output$pca_plot <- renderImage({
      outfile <- tempfile(fileext = '.png')
      list(src = file.path("pca_plot_2D.png"),
           contentType = 'image/png',
           width = 400,
           height = 300,
           alt = "This is alternate text")
    }, deleteFile = F)
    
    result <- mrMLM(fileGen = input$genofile$datapath, 
                    filePhe = file.path(dir, "results",file = "pheno.csv"),
                    fileKin = NULL, filePS = input$pcafile$datapath,
                    PopStrType = "PCA",fileCov = NULL, Genformat = "Cha",
                    method=c("mrMLM","FASTmrMLM","FASTmrEMMA","ISIS EM-BLASSO","pLARmEB"),
                    trait = 1:1, SearchRadius = 50, CriLOD = 3,SelectVariable = 50,
                    Bootstrap = FALSE, Plotformat = "jpeg",
                    dir =file.path(dir,"results"),
                    RAM = 100,DrawPlot = TRUE)
    
    #file sorting
    if(file.exists(file.path(dir,"results","1_Final result.csv"))){
      removeModal()
      
      showModal( modalDialog(
        h4(paste0("MTAs are identied and the file is present at: \nresults/1_Final result.csv")),
        footer=tagList(actionButton("ok_gwas","Proceed to download the files"))
      ))
      
      observeEvent(input$ok_gwas, {
        removeModal()
      })
      
     
      output$man <- renderImage({
        outfile <- tempfile(fileext = '.jpeg')
        list(src = file.path(dir,"results","1_Manhattan plot.jpeg"),
             contentType = 'image/jpeg',
             width = 400,
             height = 300,
             alt = "This is alternate text")
      }, deleteFile = F)
      
      output$QQ <- renderImage({
        outfile <- tempfile(fileext = '.jpeg')
        list(src = file.path(dir,"results","1_qq plot.jpeg"),
             contentType = 'image/jpeg',
             width = 400,
             height = 300,
             alt = "This is alternate text")
      }, deleteFile = F)
      
      
      #we are chaning the above to all rs ids even if they come once
      res <- read.csv(file = file.path(dir,"results","1_Final result.csv"))
      out <- res[!duplicated(res$RS.),]
      
      for (nr in c(1:nrow(out))) {
        out$Method[nr] <- toString(res$Method[res$RS.==out$RS.[nr]])
        out$QTN.effect[nr] <- mean(res$QTN.effect[res$RS.==out$RS.[nr]])
        out$r2....[nr] <- mean(res$r2....[res$RS.==out$RS.[nr]])
      }
      
      candgene_ind <- out[,c(5,6)]
      write_csv(candgene_ind, file.path(dir, "POS.csv"))
      write.csv(out[,c(2:8)],file.path(dir,"Identified_MTAs.csv"),row.names = F) 
      
      output$result_table = renderTable({
        return(out)
      })
    }else{
      output$analysis_complete = renderText({
        print("....GWAS analysis started....")
      })
    }
    
    #removeModal()
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #Hap-Phe tabs----
  observeEvent(input$runhap,{
    
    gene <- reactiveVal(NULL)
    req(input$LD)
    
    season_loc = reactive({
      req(input$hapfile)
      a <- read.csv(input$hapfile$datapath, header = T)
      return(a)
    })
    
    df <- candidate_gene(file.path(dir,"POS.csv"),input$LD, genes_file)
    df2 = df[["gene_id"]]
    write.csv(df2,"locus.csv",row.names = F)
    
    output$candtable = renderTable({
      head(df)
    })
    
    showModal( modalDialog(
      h4(paste0("Haplo-Pheno analysis for ",length(df2)," Candidate genes")),
      footer=tagList(h3("running..."))
    ))
    
    if(file.exists(file.path(dir,"locus.csv"))){
      #reactive object supplied with a forigen function
      gene <- hap_phe(gene_infile = file.path(dir,"locus.csv"),
                      pheno_file = input$hapfile$datapath,
                      select_cri = input$hl,dir)
      func_piechart(gene_infile = file.path(dir,"locus.csv"),
                    pheno_file = input$hapfile$datapath,dir)
      removeModal()
    }
    
    output$haptab = renderTable({
      head(gene)
    })
    
    showModal(modalDialog(
      h4(paste0("Haplo-Pheno analysis for is complete")),
      footer=tagList(actionButton("ok","Proceed to download the files"))
    ))
    
    observeEvent(input$ok, {
      removeModal()
    })
    
    output$canddw = downloadHandler(
      filename =  function(){
        "candidate_genes.csv"
      },
      content = function(file){
        write.csv(df,file)
      }
    )
    
    output$hapdw = downloadHandler(
      file = function(){
        "superior_haplotype.csv"
      },
      content = function(file){
        write.csv(gene,file)
      }
    )
    updateSelectInput(session,"loc_id",choices = df2)
    updateSelectInput(session,"season_pie",choices = names(season_loc()))
    
    output$piechart <- renderImage({
      outfile <- tempfile(fileext = '.png')
      
      # Return a list containing the filename
      list(src = file.path(dir,input$loc_id,input$season_pie,"hap_diversity.png"),
           contentType = 'image/png',
           width = 400,
           height = 300,
           alt = "This is alternate text")
    }, deleteFile = TRUE)
    
  })
  
  #end of server--nextline
}

shinyApp(ui,server)
