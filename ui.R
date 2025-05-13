#UI start----
ui = fluidPage(tagList(
  navbarPage(
    theme = shinythemes::shinytheme("cosmo"),
    id = "my-page",
    title = span( "HapGUI", style = "color:white"),
    tabPanel(
      "Home", icon = icon("home"),
      # Global CSS for bold h4 tags
      tags$style(HTML("
          h4 {
            font-weight: bold;
          }
        ")),
<<<<<<< HEAD
      
=======

>>>>>>> 7715484 (Cleaned history and added final files)
      fluidRow(
        column(
          width = 12,
          div(
            style = "display: flex; flex-direction: column; align-items: center; height: auto; padding-top: 50px;",
            div(
              style = "max-width: 800px; margin-bottom: 30px;",
              h3("Welcome to HapGUI", style = "text-align: center; color: #007bff;"),
              p(
<<<<<<< HEAD
                "HapGUI is a comprehensive platform designed to analyze and explore haplotype data from different germplasm sources. 
                  It allows researchers to gain insights into the genetic diversity and structure of rice populations, with a focus on 
=======
                "HapGUI is a comprehensive platform designed to analyze and explore haplotype data from different germplasm sources.
                  It allows researchers to gain insights into the genetic diversity and structure of rice populations, with a focus on
>>>>>>> 7715484 (Cleaned history and added final files)
                  both 3K-RICE and external datasets.",
                style = "text-align: justify;"
              )
            ),
            # The styled box for choices
            div(
<<<<<<< HEAD
              style = "border: 2px solid #ecf0f5; border-radius: 10px; padding: 30px; 
=======
              style = "border: 2px solid #ecf0f5; border-radius: 10px; padding: 30px;
>>>>>>> 7715484 (Cleaned history and added final files)
                         box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1); background-color: #f8f9fa; width: 400px;",
              h4("Choose Germplasm"),
              radioButtons(
                "choose_analysis",
                "Choose analysis type",
                choices = list("3K-RICE" = "option1", "External" = "option2"),
                selected = "option1",
                inline = FALSE
              ),
              uiOutput("conditional_info"),
              div(
                style = "text-align: center; margin-top: 20px;",
                actionButton("run_analysis", "Go to Analysis Tab")
              ),
            )
          )
        )
      )
    ),
<<<<<<< HEAD
    
=======

>>>>>>> 7715484 (Cleaned history and added final files)
    #---plots and BLUPs----
    tabPanel("Phenotype Analysis",fluid = T,
             
             sidebarLayout(
               sidebarPanel(
                 
                 fileInput("pheno","Upload phenotypic data",accept = ".csv"),
                 
                 #selecting input coloumns
                 selectInput("Acc","Choose genotype column",choices = NULL),
                 selectInput("Envir1","Choose 1st environment",choices = c("None" = "")), 
                 selectInput("Envir2","Choose 2nd environment",choices = c("None" = "")),
                 selectInput("Unit1","Choose 1st experiemental unit",choices = c("None" = "")), 
                 selectInput("Unit2","Choose 2nd experiemental unit",choices = c("None" = "")),
                 
                 #adding conditional UI
                 uiOutput("cond_ui_pheno"),
                 
                 actionButton("plt_bn","Generate plots"),
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("histogram plot",
                            selectInput("trait","Choose trait column",
                                        choices = NULL),
                            
                            tags$div(selectInput("Envir1_list", "Environment1", 
                                                 choices = NULL),
                                     style="display:inline-block"),
                            tags$div(selectInput("Envir2_list", "Environment2", 
                                                 choices = NULL),
                                     style="display:inline-block"),
                            tags$div(selectInput("Unit1_list", "Experimental Unit1", choices = c("None" = "")),
                                     style="display:inline-block"),
                            tags$div(selectInput("Unit2_list", "Experimental Unit2", choices = c("None" = "")),
                                     style="display:inline-block"),
                            plotOutput("hist_pl"),
                            downloadButton("hist_dw","Download_histogram")
                   ),
                   tabPanel("violin plot",
                            selectInput("trait_violin","Choose trait column",choices = c("None" = "")),
                            selectInput("violin_based","Violin plot based on:",choices = c("")),
                            plotOutput("violin_pl"),
                            downloadButton("violin_dw","Download_violin"),
                            plotOutput("other_violin"),
                            downloadButton("other_violin_dw","Download_violin")
                   ),
                   
                   tabPanel("bar plot",
                            selectInput("trait_bar","Choose trait column",choices = c("None" = "")),
                            selectInput("bar_based","Bar plot based on:",choices = c("")),
                            plotOutput("bar_pl"),
                            downloadButton("bar_dw","Download_bar"),
                            plotOutput("other_bar"),
                            downloadButton("other_bar_dw","Download_bar")
                   ),
                   
                   tabPanel("Corrplot",
                            fluidRow(
                              wellPanel(
                                h2("Multi trait correlation"),
                                selectInput("corr_envir1", "Select one Environment(1)", choices = NULL),
                                selectInput("corr_envir2", "Select one Environment(2)", choices = NULL),
                                
                                selectInput("corr_unit1", "Select one Experimental unit(1)", choices = NULL),
                                selectInput("corr_unit2", "Select one Experimental unit(2)", choices = NULL),
                                
                                p("Please select atleast 2 traits"),
                                checkboxGroupInput("corr_traits","Select traits", choices = NULL),
                                plotOutput("corr"),
                                downloadButton("cor_dw","dowload correlation plot")
                              )),
                            fluidRow(
                              wellPanel(
                                h2("Multi location & season correlation"),
                                selectInput("corr_trait", "Select a Trait", choices = NULL),
                                plotOutput("traitcorr"),
                                downloadButton("traitcor_dw","dowload correlation plot"))
                            )
                   ),
                   tabPanel("Descriptive Statsitics",
                            p("Please select the traits"),
                            checkboxGroupInput("all_traits","Select traits", choices = NULL),
                            fluidRow(h3("Overall summary")),
                            tableOutput("table"),
                            downloadButton("summary_dw","download"),
                            fluidRow(h3("Environment/Experimental Unit wise summary")),
                            tableOutput("breakup_table"),
                            downloadButton("break_summary_dw","download")
                   ),
                   
                   tabPanel("ANOVA",
                            h3("Augmented RCBD ANOVA"),
                            selectInput("anova_treatment","select treatment column",choices = NULL),
                            selectInput("anova_block","select block column",choices = NULL),
                            selectInput("anova_trait","select trait column",choices = NULL),
                            actionButton("runova","run_anova"),
                            verbatimTextOutput("anova_new"),
                   ),
                   
<<<<<<< HEAD
                   tabPanel("BLUPs",
=======
                   tabPanel("BLUP",
>>>>>>> 7715484 (Cleaned history and added final files)
                            h2("aug RCBD"),
                            tags$div(selectInput("BEnvir1_list", "Environment1", 
                                                 choices = NULL),
                                     style="display:inline-block"),
                            tags$div(selectInput("BEnvir2_list", "Environment2", 
                                                 choices = NULL),
                                     style="display:inline-block"),
                            tags$div(selectInput("BUnit1_list", "Experimental Unit1", choices = c("None" = "")),
                                     style="display:inline-block"),
                            tags$div(selectInput("BUnit2_list", "Experimental Unit2", choices = c("None" = "")),
                                     style="display:inline-block"),
                            selectInput("block2","Choose blocks column",choices = NULL),
                            checkboxGroupInput("non_trait", "Select traits", choices = NULL),
                            br(),
                            verbatimTextOutput("blupv"),
                            actionButton("genb","Generate BLUPs"),
                            br(),
                            tableOutput("blup_output"),
                            br(),
                            downloadButton("blup_dw","download data with ajusted means")
                            
                   ),
                   
                   tabPanel("Genetic Parameters",
                            selectInput("block_gav","Choose blocks column",choices = NULL),
                            selectInput("treatment_gav","select treatment column",choices = NULL),
<<<<<<< HEAD
                            checkboxGroupInput("traits_gav", "Select traits", choices = NULL),                            tableOutput("Gen_par"),
=======
                            checkboxGroupInput("traits_gav", "Select traits", choices = NULL),                            
                            # tableOutput("Gen_par"),
>>>>>>> 7715484 (Cleaned history and added final files)
                            tableOutput("Gen_par"),
                            downloadButton("Gen_dw","Download_table"),
                            #we are using index for blocks
                   )
                 )
               )
             )
<<<<<<< HEAD
      ),
   
=======
    ),
>>>>>>> 7715484 (Cleaned history and added final files)
    #Geno extract----
    tabPanel("Genofile Extract",fluid = T,
             sidebarLayout(
               sidebarPanel(
                 fileInput("geno_extract","Phenotypic file input:"),
                 radioButtons("choose_geno","Choose genotypic file",
<<<<<<< HEAD
                              choices = list("Default" = "option1", 
=======
                              choices = list("Default" = "option1",
>>>>>>> 7715484 (Cleaned history and added final files)
                                             "External" = "option2")),
                 tags$br(),
                 uiOutput("geno_cond_ui"),
                 actionButton("run_extract","Run extract")
               ),
<<<<<<< HEAD
               
=======

>>>>>>> 7715484 (Cleaned history and added final files)
               mainPanel(
                 tabsetPanel(
                   tabPanel("Principle Components Analysis",fluid = T,
                            imageOutput("pca_plot"),
                            downloadButton("plot_down","download pca plot"),
                            tableOutput("pca_table"),
                            downloadButton("table_down","download pca table"),
                   ),
                   tabPanel("Geno",fluid = T,
                            tableOutput("geno_table"),
                            downloadButton("geno","download genotypic data"),
                   ),
                   tabPanel("Annotation",
                            fluidRow(
                              wellPanel(
                                h2("Annotated variants"),
                                tableOutput("anno_table"),
                                downloadButton("table_down","download annotation results")
                              )
                            )
                   ),
                 )
               )
             )
    ),
<<<<<<< HEAD
    
=======

>>>>>>> 7715484 (Cleaned history and added final files)
    # GWAS tab ----------------------------------------------------------------
    tabPanel("GWAS",fluid = T,
             sidebarLayout(
               sidebarPanel(
                 fileInput("gwasfile","Phenotypic file input:"),
                 numericInput("pca_obs", "Principal Componenets", 3, min = 0, max = 10),
                 checkboxGroupInput("gwas_methods", "Choose GWAS methods:",
                                    choices = list("mrMLM" = "mrMLM", "FASTmrMLM" = "FASTmrMLM",
                                                   "FASTmrEMMA" = "FASTmrEMMA", "ISIS EM-BLASSO" = "ISIS EM-BLASSO",
                                                   "pLARmEB" = "pLARmEB"),
                                    selected = "mrMLM"), # Default selection
                 #actionButton("select_all_methods", "Select All"),
                 actionButton("runa","Run GWAS")
               ),
               mainPanel(
                 tabsetPanel(
                   #------Association
                   tabPanel("Association",
                            plotOutput("man"),
                            plotOutput("QQ")
                   ),
<<<<<<< HEAD
                   
=======

>>>>>>> 7715484 (Cleaned history and added final files)
                   #----------MTAs
                   tabPanel("MTAs",fluid = T,
                            tableOutput("result_table"),
                            downloadButton("downloadData","Download GWAS results")
                   )
                 )
               )
             )
    ),
<<<<<<< HEAD
    
=======

>>>>>>> 7715484 (Cleaned history and added final files)
    # Et-GWAS tab ----
    tabPanel("Et-GWAS",fluid = T,
             sidebarLayout(
               sidebarPanel(
                 selectInput(inputId = "Bulk_size",
                             label = "Bulk size:",
                             choices = c("10", "15", "20")),
                 
                 textInput("Trait", "Trait" , "Ex: SPY "),
                 verbatimTextOutput("value"),
                 fileInput("etgwas_pheno", "Phenotype"),
                 actionButton("run_etgwas", "Run"),
                 downloadButton("downloadEt_Data", label = "Download"),
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Phenotypic distribution",
                            plotOutput("plot1", width = "100%",
                                       height = "300px"),
                            plotOutput("plot2", width = "100%", height = "300px")),
                   tabPanel("Association",
                            plotOutput("plot3", width = "100%",
                                       height = "300px")))
               )
             )
    ),
    # hap-phe tab -------------------------------------------------------------
    tabPanel("Haplo-Pheno",fluid = T,
             sidebarLayout(
               sidebarPanel(
                 fileInput("hapfile","Phenotypic file input:"),
                 fileInput("posfile","MTA file input:"),
<<<<<<< HEAD
                 numericInput("LD", "LD region", 25000,min = 1000, max = 250000),
                 radioButtons("select_cri", "Choose High/Low Value Haplotypes", 
                              choices = c("High", "Low")),
                 actionButton("runhap","Run Haplo-Pheno")
                 
=======
                 # fileInput("gff_file", "Upload GFF File", accept = ".gff3"),
                 numericInput("LD", "LD region", 25000,min = 1000, max = 250000),
                 radioButtons("select_cri", "Choose High/Low Value Haplotypes",
                              choices = c("High", "Low")),
                 actionButton("runhap","Run Haplo-Pheno")

>>>>>>> 7715484 (Cleaned history and added final files)
               ),
               mainPanel(
                 tabsetPanel(
                   #------Candidate genes
                   tabPanel("Candidate Genes",
                            tableOutput("candtable"),
                            br(),
                            downloadButton("canddw")
                   ),
                   tabPanel("Superior haplotypes",
                            tableOutput("haptab"),
                            br(),
                            downloadButton("hapdw","Download haplotype table"),
                            br(),
                            htmlOutput('text1'),
                            br(),
<<<<<<< HEAD
                            fluidRow(
                              column(6,selectInput("loc_id","Choose Locus ID",
                                                   choices = NULL)),
                              column(6,selectInput("season_pie","Choose the season",
                                                   choices = NULL))),
                            fluidRow(
                              column(6,imageOutput("piechart")),
                              column(6,tableOutput("donortab"))),
=======
                            # fluidRow(
                            #   column(6,selectInput("loc_id","Choose Locus ID",
                            #                        choices = NULL)),
                            #   column(6,selectInput("season_pie","Choose the season",
                            #                        choices = NULL))),
                            fluidRow(
                            column(6, selectInput("loc_id", "Choose Locus ID", choices = NULL)),
                            # column(6,selectInput("season_pie","Choose the season",choices = NULL))
                            ),
                            fluidRow(
                              column(6,imageOutput("piechart")),
                              column(6,tableOutput("donortab"))
                              ),
>>>>>>> 7715484 (Cleaned history and added final files)
                            downloadButton("downloadhapData")
                   ),
                 )
               )
             )
    )
<<<<<<< HEAD
    
    
  )
)
)
=======


  )
))
>>>>>>> 7715484 (Cleaned history and added final files)
