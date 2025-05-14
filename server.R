server <- function(input, output, session) {

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
  mysummary <- function(x,na.rm=T){
    res <- list(Minimum=min(x, na.rm=na.rm),
                Maximum=max(x, na.rm=na.rm),
                Mean=mean(x, na.rm=na.rm),
                StdDev=sd(x,na.rm=na.rm),
                Sknewness = skewness(x),
                Kurtosis = kurtosis(x))
    
    res$Range <- paste0(res$Minimum,"-",res$Maximum)
    res
  }
  # Conditional information for each choice
  output$conditional_info <- renderUI({
    if (input$choose_analysis == "option1") {
      tagList(
        p("The 3K-RICE dataset includes over 3,000 rice genomes, providing a rich resource for studying genetic variation.",
          style = "text-align: justify; color: #333;")
      )
      
    } else if (input$choose_analysis == "option2") {
      tagList(
        p("Upload your external germplasm data to analyze its haplotypes and compare it with existing datasets.",
          style = "text-align: justify; color: #333;")
      )
    } else {
      NULL
    }
  })
  
  observeEvent(input$run_analysis, {
    print("Button clicked!")  # Debugging
    updateNavbarPage(session, "my-page", selected = "Phenotype Analysis")
  })
  
  #update Pheno tab
  output$cond_ui_pheno <- renderUI({
    req(input$choose_analysis)  # Ensure input$choose_analysis is available
    
    if (input$choose_analysis == "option1") {
      tagList(
        selectInput("type", "Select genotype ID", choices = c("Accession", "IRIS_ID", "IRGC_NO"))
      )
    } else {
      tagList(
        selectInput("category","Choose category column",
                    choices = c("None" = "", names(trait()))),
        tags$br()
      )
    }
  })
  
  pop =  reactive({
    req(input$choose_analysis)
    if(input$choose_analysis == "option1"){
      pop <- read.delim("Edited all population.txt", header = T)
      names(pop) <- c("Designation","Name","Subpopulation","COUNTRY","IRGC")
      return(pop)
    }
    return(data.frame())
  })
  
  #Taking phenotypic file as reactive input
  trait = reactive({
    req(input$pheno)
    a <- read.csv(input$pheno$datapath, header = T)
    return(a)
  })
  
  observeEvent(trait(),{
    req(trait())
    
    updateSelectInput(session,"Acc", choices = names(trait()))
    updateSelectInput(session,"Envir1",choices = c("None" = "", names(trait())))
    updateSelectInput(session,"Envir2",choices = c("None" = "", names(trait())))
    updateSelectInput(session,"Unit1",choices = c("None" = "", names(trait())))
    updateSelectInput(session,"Unit2",choices = c("None" = "", names(trait())))
    #updateCheckboxGroupInput(session,"traits",choices = names(trait()))
  })
  
  observeEvent(input$plt_bn,{
    req(trait())
    input_colnames = names(trait())
    envir1 = unique(trait()[[input$Envir1]])
    envir2 = unique(trait()[[input$Envir2]])
    
    unit1 = unique(trait()[[input$Unit1]])
    unit2 = unique(trait()[[input$Unit2]])
    
    #update histogram tab
    updateSelectInput(session, "trait", 
                      choices = c("None", names(trait())),selected = "None")
    updateSelectInput(session,"Envir1_list",choices = append("None",envir1))
    updateSelectInput(session,"Envir2_list",choices = append("None",envir2))
    
    updateSelectInput(session,"Unit1_list",choices = append("None",unit1))
    updateSelectInput(session,"Unit2_list",choices = append("None",unit2))
    
    #update violin tab
    updateSelectInput(session,"trait_violin",
                      choices = c("None", names(trait())),selected = "None")
    updateSelectInput(session,"violin_based",
                      choices = c(input$Envir1,input$Envir2,input$Unit1,input$Unit2))
    
    #update bar tab
    updateSelectInput(session,"trait_bar",
                      choices = c("None", names(trait())),selected = "None")
    updateSelectInput(session,"bar_based",
                      choices = c(input$Envir1,input$Envir2,input$Unit1,input$Unit2))
    
    #corr
    updateSelectInput(session,"corr_envir1",choices = append("None",envir1))
    updateSelectInput(session,"corr_envir2",choices = append("None",envir2))
    
    updateSelectInput(session,"corr_unit1",choices = append("None",unit1))
    updateSelectInput(session,"corr_unit2",choices = append("None",unit1))
    
    updateCheckboxGroupInput(session,"corr_traits",choices = names(trait()))
    
    updateSelectInput(session,"corr_trait",choices = append("None",names(trait())))
    
    #update descriptive tab
    updateCheckboxGroupInput(session,"all_traits",choices = names(trait()))
    
    #new anova update
    updateSelectInput(session,"anova_treatment",choices = names(trait()))
    updateSelectInput(session,"anova_block",choices = names(trait()))
    updateSelectInput(session,"anova_trait",choices = names(trait()))
    
    #update blup tab
    updateSelectInput(session,"BEnvir1_list",choices = append("None",envir1))
    updateSelectInput(session,"BEnvir2_list",choices = append("None",envir2))
    
    updateSelectInput(session,"BUnit1_list",choices = append("None",unit1))
    updateSelectInput(session,"BUnit2_list",choices = append("None",unit2))
    
    updateSelectInput(session,"block2",choices = names(trait()))
    updateCheckboxGroupInput(session,"non_trait",choices = names(trait()))
    
    #update GAV tab
    updateSelectInput(session,"block_gav",choices = names(trait()))
    updateSelectInput(session,"treatment_gav",choices = names(trait()))
    updateCheckboxGroupInput(session,"traits_gav",choices = names(trait()))
    
  })
  
  data <- reactive({
    req(input$pheno, trait(),input$plt_bn)  
    # Switch logic based on condition
    if (input$choose_analysis == "option1") {
      trait_d = trait()
      names(trait_d)[which(names(trait_d)== input$Acc)] <- c ("Name")
      
      if (!is.null(input$trait)){
        trait_d$phenotype <- trait_d[[input$trait]]
      }
      
      if (!is.null(input$trait_violin)) {
        trait_d$phenotype_violin <- trait_d[[input$trait_violin]]
      }
      
      if (!is.null(input$trait_bar)) {
        trait_d$phenotype_bar <- trait_d[[input$trait_bar]]
      }
      
      ycol <- "Name"
      if(input$type == "IRIS_ID"){
        ycol <- "Designation"
        trait_d$Name <- gsub("IRIS 313","IRIS_313",trait_d$Name)
      }else if(input$type == "IRGC_NO"){
        ycol <- "IRGC"
      }
      aa <- merge(trait_d, pop(), by.x = "Name", by.y = ycol, all = F)
      subpop <- c("indx", "ind2", "ind1B", "ind1A", "ind3", "aus", "japx", "temp", "trop", "subtrop", "admix", "aro")
      subcol <- c("ind", "ind", "ind", "ind", "ind", "aus", "jap", "jap", "jap", "jap", "admix", "aro")
      aa$sub <- aa$Subpopulation
      for (i in c(1:12)) {
        aa$sub[aa$Subpopulation == subpop[i]] <- subcol[i]
      }
      
    } else {
      trait_d = trait()
      names(trait_d)[which(names(trait_d)== input$Acc)] <- c ("Name")
      
      if (!is.null(input$trait)){
        trait_d$phenotype <- trait_d[[input$trait]]
      }
      
      if (!is.null(input$trait_violin)) {
        trait_d$phenotype_violin <- trait_d[[input$trait_violin]]
      }
      
      if (!is.null(input$trait_bar)) {
        trait_d$phenotype_bar <- trait_d[[input$trait_bar]]
      }
      
      if(!is.null(input$category) & input$category != "None"){
        names(trait_d)[which(names(trait_d)== input$category)] <- c ("subpop")
      }
      aa = trait_d
    }
    return(aa)
  })
  
  
  #---Histogram
  histo =  reactive({
    req(data())
    pheno_d2 = data()
    
    if(!is.null(input$Envir1_list) & input$Envir1_list != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Envir1]]) == input$Envir1_list, ]
    }
    
    if(!is.null(input$Envir2_list) & input$Envir2_list != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Envir2]]) == input$Envir2_list, ]
    }
    
    if(!is.null(input$Unit1_list) & input$Unit1_list != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Unit1]]) == input$Unit1_list, ]
    }
    
    if(!is.null(input$Unit2_list) & input$Unit2_list != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Unit2]]) == input$Unit2_list, ]
    }
    
    p1 <- ggplot(data= pheno_d2, aes(phenotype)) + theme_bw() +
      theme(axis.line = element_line(size=1, colour = "black"),
            panel.grid = element_blank(), panel.border = element_blank(),
            panel.background = element_blank(), panel.spacing = unit(0.3, "lines"),
            axis.text = element_text(colour="black", size = 12,family = "Arial"),
            #axis.text.y=element_text(colour="black", size = 12),
            axis.title=element_text(size=12,face = "bold",color ="black"),
            strip.text.x = element_text(size = 12, face="bold"))+
      geom_histogram(data=pheno_d2,color="black", fill = "skyblue",
                     alpha = 0.5)  + xlab(input$trait) + ylab("Count")
    return(p1)
  })
  
  output$hist_pl = renderPlot({
    req(data(),input$trait)
    histo()
  })
  
  output$hist_dw = downloadHandler(
    filename = function(){
      paste0(input$trait," histogram_plot.png")
    },
    content = function(file){
      ggsave(file,histo(), width = 8, height = 6,dpi = 600)
    }
  )
  
  #---Violin
  violin = reactive({
    req(data())
    p2 <- ggplot(data(), aes(x = sub,y=phenotype_violin, fill=sub)) +
      geom_violin()+ theme_cus +
      ylab(input$trait_violin) + geom_boxplot(width=0.1) +
      scale_fill_brewer()+xlab(" ")
    return(p2)
  })
  
  output$violin_pl = renderPlot({
    req(data())
    violin()
  })
  
  output$violin_dw = downloadHandler(
    filename = function(){
      paste0(input$trait_violin," violin_plot.png")
    },
    content = function(file){
      ggsave(file,violin(), width = 8, height = 6,dpi = 600)
    }
  )
  
  violin_others <- reactive({
    req(data(),input$violin_based)
    pheno_d2 = data()
    
    p2 <- ggplot(pheno_d2, aes(x = as.factor(pheno_d2[[input$violin_based]]),
                               y = pheno_d2$phenotype_violin,
                               fill = factor(pheno_d2[[input$violin_based]]))) + 
      geom_violin() + labs(x = input$violin_based,
                           y = input$trait_violin)+
      scale_fill_discrete(name = input$violin_based)
    
    return(p2)
  })
  
  output$other_violin = renderPlot({
    req(input$plt_bn)
    violin_others()
  })
  
  output$other_violin_dw = downloadHandler(
    filename = function(){
      paste0(input$trait_violin,"_",input$violin_based," violin_plot.png")
    },
    content = function(file){
      ggsave(file,violin_others(), width = 8, height = 6,dpi = 600)
    }
  )
  
  #barplot
  bar = reactive({
    req(data())
    
    df <- data_summary(data(), varname="phenotype_bar",
                       groupnames="sub")
    
    colnames(df) <- c("sub","phenotype","sd")
    print("summary table")
    
    p3 <- ggplot(df, aes(x = sub,y=phenotype, fill=sub)) +
      geom_bar(stat="identity", color="black",
               position=position_dodge()) +
      geom_errorbar(aes(ymin=phenotype-sd, ymax=phenotype+sd), width=.2,
                    position=position_dodge(.9)) +
      theme_cus + ylab(input$trait_bar)+
      scale_fill_brewer()+xlab(" ")
    
    return(p3)
  })
  
  output$bar_pl = renderPlot({
    req(input$plt_bn)
    bar()
  })
  
  output$bar_dw = downloadHandler(
    filename = function(){
      paste0(input$trait_bar,"_barplot.png")
    },
    content = function(file){
      ggsave(file,bar(), width = 8, height = 6, dpi = 600)
    }
  )
  
  bar_others = reactive({
    req(data(),input$bar_based)
    pheno_d2 = data()
    df <- data_summary(data(), varname="phenotype_bar",
                       groupnames=input$bar_based)
    
    colnames(df) =c(input$bar_based,"phenotype","sd")
    
    p3 <- ggplot(df, aes(x = df[[input$bar_based]],y=phenotype, 
                         fill = factor(df[[input$bar_based]]))) +
      geom_bar(stat="identity", color="black",
               position=position_dodge()) +
      geom_errorbar(aes(ymin=phenotype-sd, ymax=phenotype+sd), width=.2,
                    position=position_dodge(.9)) + theme_cus + 
      ylab(input$trait_bar)+ scale_fill_discrete(name = input$bar_based) +
      xlab("")
    
    return(p3)
  })
  
  output$other_bar <- renderPlot({
    req(input$plt_bn, input$bar_based)
    bar_others()
  })
  
  output$other_bar_dw = downloadHandler(
    filename = function(){
      paste0(input$trait_bar,"_",input$bar_based,"_barplot.png")
    },
    content = function(file){
      ggsave(file,bar_others(), width = 8, height = 6, dpi = 600)
    }
  )
  
  #correlation----
  corr = reactive({
    req(trait())
    pheno_d2 = trait()
    
    if(!is.null(input$corr_envir1) & input$corr_envir1 != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Envir1]]) == input$corr_envir1,]
    }
    if(!is.null(input$corr_envir2) & input$corr_envir2 != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Envir2]]) == input$corr_envir2,]
    }
    if(!is.null(input$corr_unit1) & input$corr_unit1 != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Unit2]]) == input$corr_unit1,]
    }
    if(!is.null(input$corr_unit2) & input$corr_unit2 != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Unit2]]) == input$corr_unit2,]
    }
    
    dd <- pheno_d2[,input$corr_traits]
    dd_clean <- na.omit(dd)
    
    png(file.path(dir,"multitrait_corr.png"),
        width = 6, height = 6, units = "in", res = 300)
    
    pairs.panels(dd_clean, smooth = TRUE, scale = FALSE,
                 density = TRUE, ellipses = TRUE, method = "pearson",
                 pch = 21, lm = FALSE, cor = TRUE, jiggle = FALSE,
                 factor = 2, hist.col = 3, stars = TRUE, ci = TRUE)
    dev.off()
  })
  
  output$corr <- renderImage({
    outfile <- tempfile(fileext = '.png')
    corr()
    # Return a list containing the filename
    list(src = file.path(dir,"multitrait_corr.png"),
         contentType = 'image/png',
         width = 400,
         height = 400,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  output$cor_dw = downloadHandler(
    filename <- function() {
      "multitrait_corr.png"
    },
    content <- function(file) {
      file.copy(file.path(dir,"multitrait_corr.png"), file)
    },
    contentType = "image/png"
  )
  
  corr2 = reactive({
    selected_cols <- c(input$Envir1, input$Envir2, input$Unit1, 
                       input$Unit2)
    selected_cols <- selected_cols[selected_cols != "None" & selected_cols != ""]
    
    aa <- trait()
    aa <- aa[,c(selected_cols,input$Acc,input$corr_trait)]
    print(head(aa))
    aa <- aa %>%
      mutate(Combined = unite(aa, "Combined", 
                              names(aa)[1:(ncol(aa) - 1)], sep = "_")$Combined)
    
    aa <- aa[!duplicated(aa$Combined),]
    aa <- aa[,-c(ncol(aa))]
    print(head(aa))
    reshaped_data <- aa %>%
      unite("Combined", names(aa)[1:(ncol(aa) - 2)], sep = "_")
    
    print(head(reshaped_data))
    
    df <- reshaped_data %>% spread(key = Combined, value = input$corr_trait)
    head(df)
    mat <- df[,-c(1)]
    mat <- na.omit(mat)
    
    png(file.path(dir,paste0(input$corr_trait,"_corr.png")),
        width = 6, height = 6, units = "in", res = 300)
    
    pairs.panels(mat, smooth = TRUE, scale = FALSE,
                 density = TRUE, ellipses = TRUE, method = "pearson",
                 pch = 21, lm = FALSE, cor = TRUE, jiggle = FALSE,
                 factor = 2, hist.col = 3, stars = TRUE, ci = TRUE)
    dev.off()
  })
  
  output$traitcorr <- renderImage({
    outfile <- tempfile(fileext = '.png')
    corr2()
    # Return a list containing the filename
    list(src = file.path(dir,paste0(input$corr_trait,
                                    "_corr.png")),
         contentType = 'image/png',
         width = 400,
         height = 400,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  output$traitcor_dw = downloadHandler(
    filename <- function() {
      paste0(input$corr_trait,"_corr.png")
    },
    content <- function(file) {
      file.copy(file.path(dir,paste0(input$trait,
                                     "_corr.png")), file)
    },
    contentType = "image/png"
  )
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  summ = reactive({
    req(trait(),input$Acc)
    column_data <- trait()[,c(input$Acc,input$all_traits)]
    
    colnames(column_data)[1] <- "Name"
    
    print(head(column_data))
    melted_data <- setDT(melt(setDT(column_data),
                              id.vars = c("Name")))
    print(head(melted_data))
    summary_stats <- melted_data[,mysummary(value),
                                 by=.(variable)]
    print(summary_stats)
    setnames(summary_stats, "variable", "Trait")
    return(summary_stats)
  })
  
  break_summ = reactive({
    req(trait())
    selected_cols <- c(input$Envir1, input$Envir2, input$Unit1, 
                       input$Unit2)
    
    selected_cols <- selected_cols[selected_cols != "None" & selected_cols != ""]
    
    column_data <- trait()[,c(selected_cols,input$all_traits)]
    melted_data <- setDT(melt(setDT(column_data),
                              id.vars=selected_cols))
    summary_stats <- melted_data[, mysummary(value), 
                                 by=c("variable", selected_cols)]
    
    setnames(summary_stats, "variable", "Trait")
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
      "Envir_Unit_summary.csv"
    },
    content = function(file){
      write.csv(break_summ(),file)
    }
  )
  
  #ANOVA: augRCBD
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
  
  #BLUPs---
  blup_out = reactive({
    req(input$genb, trait(), data())
    pheno_d2 = trait()
    
    if(!is.null(input$BEnvir1_list) & input$BEnvir1_list != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Envir1]]) == input$BEnvir1_list, ]
    }
    
    if(!is.null(input$BEnvir2_list) & input$BEnvir2_list != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Envir2]]) == input$BEnvir2_list, ]
    }
    
    if(!is.null(input$BUnit1_list) & input$BUnit1_list != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Unit1]]) == input$BUnit1_list, ]
    }
    
    if(!is.null(input$BUnit2_list) & input$BUnit2_list != "None"){
      pheno_d2 <- pheno_d2[as.character(pheno_d2[[input$Unit2]]) == input$BUnit2_list, ]
    }
    
    d <- pheno_d2
    names(d)[which(names(d)== input$Acc)] <- c ("Name")
    d = na.omit(d)
    #making blocks as factors
    d[[input$block2]] = as.factor(d[[input$block2]])
    
    #prepare terms for BLUP calc
    acc_group <- d[["Name"]]
    block = d[[input$block2]]
    
    l = list()
    for (i in input$non_trait){
      model<-lmer(d[[i]]~(1|acc_group)+(1|block),data = d)
      blup_val <- as.data.frame(coef(model)$acc_group)
      names(blup_val)[which(names(blup_val)=="(Intercept)")] = paste0(i,"(adjusted)")
      l[[paste0(i,"(adjusted)")]] = blup_val
    }
    
    k = bind_cols(l)
    k <- k %>% dplyr::select(ends_with("(adjusted)"))
    k$"Genotype" = row.names(k)
    print("the blup first data")
    print(head(k))
    
    if (input$choose_analysis == "option1") {
      names(d)[which(names(d)== input$Acc)] <- c ("Name")
      
      ycol <- "Name"
      if(input$type == "IRIS_ID"){
        ycol <- "Designation"
        d$Name <- gsub("IRIS 313","IRIS_313",d$Name)
      }else if(input$type == "IRGC_NO"){
        ycol <- "IRGC"
      }
      
      aas <- merge(d, pop(), by.x = "Name", by.y = ycol, all = F)
      subpop <- c("indx", "ind2", "ind1B", "ind1A", "ind3", "aus", "japx", "temp", "trop", "subtrop", "admix", "aro")
      subcol <- c("ind", "ind", "ind", "ind", "ind", "aus", "jap", "jap", "jap", "jap", "admix", "aro")
      aas$sub <- aas$Subpopulation
      for (i in c(1:12)) {
        aas$sub[aas$Subpopulation == subpop[i]] <- subcol[i]
      }
      
    } else {
      d = trait()
      names(d)[which(names(d)== input$Acc)] <- c ("Name")
      
      if(!is.null(input$category) & input$category != "None"){
        names(d)[which(names(d)== input$category)] <- c ("subpop")
      }
      aas = d
    }                             
    
    tot_col <- ncol(d)+1
    aas <- aas[,c(tot_col,1:ncol(d),(tot_col+1):ncol(aas))]
    
    final =  merge(aas,k,by.x = "Name",by.y = "Genotype")
    
    if (input$choose_analysis == "option1") {
      ycol <- "Name"
      if(input$type == "IRIS_ID"){
        ycol <- "Designation"
      }else if(input$type == "IRGC_NO"){
        ycol <- "IRGC"
      }
      if(ycol=="Designation"){
        names(final)[which(names(final)== "Name")] <- c ("X.Phenotype.")
      } else{
        names(final)[which(names(final)== "Designation")] <- c ("X.Phenotype.")
      }
    } else{
      names(final)[which(names(final)== "Name")] <- c ("X.Phenotype.")
    }
    return(final)
  })
  
  output$blup_output = renderTable({
    final <- blup_out()
    return(head(final))
  })
  
  output$blup_dw = downloadHandler(
    filename = function(){
      "blups.csv"
    },
    content = function(file){
      write.csv(blup_out(),file)
    }
  )
  
  #creating Genetic parameter reactive for downloads---
  GAV = reactive({
    req(input$Acc, trait())
    
    d = trait()
    d <- na.omit(d)
    
    selected_cols <- c(input$Envir1, input$Envir2, input$Unit1, 
                       input$Unit2)
    selected_cols <- selected_cols[selected_cols != "None" & selected_cols != ""]
    d[selected_cols] <- lapply(d[selected_cols], as.factor)
    
    #experiment using names(df)
    out <- augmentedRCBD.bulk(d, block = input$block_gav,
                              treatment = input$treatment_gav, 
                              traits = input$traits_gav,
                              freqdist = TRUE, gva = TRUE,
                              console = TRUE)
    return(out[["Genetic variability analysis"]])
  })
  
  output$Gen_par = renderTable({
    print(GAV())
  })
  #GAV downloads
  
  output$Gen_dw = downloadHandler(
    filename = function(){
      "Genetic Parameters.csv"
    },
    content = function(file){
      write.csv(GAV(),file)
    })
  

  #Genomic extract ----
  pheno_ingeno = reactive({
    req(input$geno_extract$datapath)
    a <- read.csv(input$geno_extract$datapath, header = T)
    return(a)
  })

  output$geno_cond_ui <- renderUI({
    req(pheno_ingeno())
    if (input$choose_geno == "option2") {
      tagList(
        p("Please upload the genotypic files in VCF format, below:",
          style = "color: blue;"),
        fileInput("geno_vcf", "Upload VCF File", accept = c(".vcf", ".vcf.gz")),
        fileInput("gff_file", "Upload GFF File", accept = ".gff3"),
        selectInput("category_geno","Choose category column",
                    choices = c("None" = "", names(pheno_ingeno())))
        # selectInput("category", "Choose grouping column", choices = NULL, selected = NULL)
      )}
  })

  ##-------- change
  # observe({
  #   req(pheno_ingeno())
  #   phe <- pheno_ingeno()
  #   extra_cols <- setdiff(colnames(phe), c("ID", "trait"))
  #   
  #   if (length(extra_cols) > 0) {
  #     updateSelectInput(session, "category", choices = extra_cols)
  #   } else {
  #     updateSelectInput(session, "category", choices = NULL)
  #   }
  # })
  #-------------
  
  data_ingeno <- reactive({
    req(input$geno_extract, pheno_ingeno(),input$run_extract)

    phe = pheno_ingeno()

    out1 <- phe[, c(1, 1)]
    write.table(out1, file = "id.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

    # Switch logic based on condition
    if (input$choose_geno == "option1") {
      colnames(phe) <- c("ID", "trait")
      phe$ID <- gsub(pattern = "IRIS ", replacement = "IRIS_", phe$ID)
      
      pop <- read.delim(subpop, header = TRUE)
      colnames(pop) <- c("ID", "subpop")
      print("Merging phenotype with subpop...")
      comb_phe <- merge(phe, pop, by = "ID", all.x = TRUE)
      # comb_phe <- phe   ####### change
    } else if (input$choose_geno == "option2") {
      
      print("Column names before renaming:")
      print(colnames(phe))
      colnames(phe)[1:2] <- c("ID", "trait")
      
      if (!is.null(input$category) && input$category != "None") {
        print("Selected category from dropdown:")
        print(input$category)
        
        if (input$category %in% colnames(phe)) {
          colnames(phe)[which(colnames(phe) == input$category)] <- "subpop"
        } else {
          warning("Selected category column does not exist in phenotype file.")
          phe$subpop <- "Unknown"
        }
      } else {
        warning("No category selected. Using default 'Unknown'.")
        phe$subpop <- "Unknown"
      }
      comb_phe <- phe
      comb_phe$subpop <- factor(comb_phe$subpop)
      
      print("Final column names:")
      print(colnames(phe))
    }
    
    # })
    
  #   # Check if 'sub' column exists
    if (!"subpop" %in% colnames(comb_phe)) {
      stop("Error: 'subpop' column is missing in comb_p")
    }
    #comb_phe should have 3 columns: accessions, phenotypic value, sub
    return(comb_phe)
  })

  observeEvent(input$run_extract, {
    showModal(modalDialog(h4("Genotypic data is being extracted")))
    # req(data_ingeno())
    phefile = data_ingeno()
    
    # phefile = comb_phe  # ERROR

    # gff_path <- NULL

    if (input$choose_geno == "option1") {
      genofile <- "marker_main"
      gff_path = paste0(dir,"/Oryza_sativa.IRGSP-1.0.60.gff3")
      genome_name <- sub(".*/([^./]+)\\..*", "\\1", gff_path)
      system(paste0(ip_dir, "/plink2 --bfile ", genofile,
                    " --keep id.txt --export vcf --out marker"))
      system(paste0(ip_dir, "/plink2 --bfile ", genofile, 
                    " --keep id.txt --make-bed --out marker_etgwas"))

    } else if (input$choose_geno == "option2" && !is.null(input$geno_vcf)) {
      #vcf_path <- input$vcf_file$datapath
      print("pringting GFF...")
      req(input$gff_file$datapath)
      gff_path <- input$gff_file$datapath
      print(gff_path)
      print(head(gff_path))
      genome_name <- sub(".*/([^./]+)\\..*", "\\1", gff_path)
      vcf <- read.vcfR(input$geno_vcf$datapath)
      print(head(vcf))
      sample_ids <- phefile$ID
      vcf_samples <- colnames(as.data.frame(vcf@gt))
      matching_samples <- intersect(sample_ids, vcf_samples)
      if (length(matching_samples) == 0) {
        stop("No matching samples found between VCF and phenotype file.")
      }

      vcf_fix <- as.data.table(vcf@fix)
      vcf_fix$ID <- ifelse(vcf_fix$ID == ".",
                           paste0(vcf_fix$CHROM, "_", vcf_fix$POS), vcf_fix$ID)
      vcf@fix <- as.matrix(vcf_fix)
      vcf_subset <- vcf[, c("FORMAT", matching_samples)]
      write.vcf(vcf_subset, file = "marker.vcf")
      # system(paste0(ip_dir, "/plink2 --vcf marker.vcf --export vcf --out marker"))
      system(paste0(ip_dir, "/plink2 --vcf marker.vcf --export vcf --out marker"))
      system(paste0(ip_dir, "/plink2 --vcf marker.vcf --make-bed --out marker_etgwas --double-id"))
    }

    choice_geno = input$choose_geno
    marker <- qc_linux(choice_geno,phefile = phefile,ip_dir = ip_dir, theme_cus)
    
    # genome_name <- paste0(sub("\\..*", "", gff_path))

    annotation_result <- vcf_annotation("marker.vcf", gff_path,
                                        ann_path, genome_name)

    removeModal()

    output$pca_plot <- renderImage({
      list(
        src = file.path(dir, "pca_plot_2D.png"),
        contentType = 'image/png',
        width = 500, height = 400,
        alt = "PCA plot"
      )
    }, deleteFile = FALSE)

    output$plot_down <- downloadHandler(
      filename <- function() {
        paste("pca", "png", sep=".")
      },
      content <- function(file) {
        file.copy("pca_plot_2D.png", file)
      },
      contentType = "image/png"
    )

    output$pca_table <- renderTable({
      req(file.exists(file.path(dir,"pca.csv")))
      read.csv(file.path(dir,"pca.csv")) %>% head()
    })

    output$geno_table <- renderTable({
      head(marker)
    })

    output$geno <- downloadHandler(
      filename = function() { "marker.csv" },
      content = function(file) { write.csv(marker, file, row.names = FALSE) }
    )

    output$table_down <- downloadHandler(
      filename = function() { "pca.csv" },
      content = function(file) { file.copy("pca.csv", file) },
      contentType = "application/csv"
    )

    output$anno_table <- renderTable({
      req(file.exists(file.path(dir,"annotated_variants.csv")))
      head(read.csv(file.path(dir,"annotated_variants.csv")), 10)
    })

    output$table_down <- downloadHandler(
      filename = function() { "annotated_variants.csv" },
      content = function(file) {
        file.copy("annotated_variants.csv", file)
      }
    )

  })

  #GWAS tabs----
  output$analysis_started = renderPrint({
    req(input$runa)
    cat("analysis started","\n")
  })

  plotData <- reactiveVal(NULL)

  observeEvent(input$runa, {
    req(input$gwasfile)
    methods_selected <- input$gwas_methods
    updateCheckboxGroupInput(session, "gwas_methods",
                             selected = c("mrMLM", "FASTmrMLM", "FASTmrEMMA",
                                          "ISIS EM-BLASSO", "pLARmEB"))

    phe = reactive({
      phe <- read.csv(file = input$gwasfile$datapath,header = T)
      # if (input$choose_geno == "option1") {
      #   phe <- phe
      # } else if (input$choose_geno == "option2") {
        phe <- phe[, 1:2]
      # } else {
      #   stop("Error: No GFF file provided!")
      # }
      pnam <- c("X.Phenotype.",names(phe)[2]) #Change trait
      colnames(phe) <- pnam
      # phe$"X.Phenotype." <- gsub(pattern = "IRIS ",replacement = "IRIS_",
      #                            phe$"X.Phenotype.")
      return(phe)
    })

    showModal( modalDialog(
      h4(paste0("GWAS for ",names(phe())[2],":")),
      footer=tagList(h3("running..."))
    ))

    if (!dir.exists(file.path(dir, "GWAS_results"))) {
      dir.create(file.path(dir, "GWAS_results"))
    }
    write.csv(phe() ,file.path(dir,"GWAS_results",file = "pheno.csv"),row.names = F)
    pca <- read.csv(file.path(dir,"pca.csv"), header = T)
    cn <- colnames(pca)
    pca1 <- insertRows(pca, 1 , new = NA)
    pca1[1,] <- cn
    pca1[1,1] <- "<ID>"
    # pca1 <- pca[,c(1:(as.numeric(input$pca_obs)+1))]
    colnames(pca1) <- "<PCA>"
    colnames(pca1)[2:ncol(pca1)] <- ""
    write.csv(pca1, file = "pca_temp.csv", row.names = FALSE)
    # Run mrMLM with the selected methods
    result <- mrMLM(fileGen = file.path(dir,"marker.csv"),
                    filePhe = file.path(dir,"GWAS_results",file = "pheno.csv"),
                    fileKin = NULL, filePS = file.path(dir,"pca_temp.csv"),
                    PopStrType = "PCA", fileCov = NULL, Genformat = "Cha",
                    method = methods_selected, # dynamically pass selected methods here
                    trait = 1:1, SearchRadius = 50, CriLOD = 3,
                    SelectVariable = 50, Bootstrap = FALSE,
                    Plotformat = "jpeg", dir = file.path(dir,"GWAS_results"),
                    RAM = 100, DrawPlot = TRUE)

    # ...

    #file sorting
    if(file.exists(file.path(dir,"GWAS_results","1_Final result.csv"))){
      removeModal()

      showModal( modalDialog(
        h4(paste0("MTAs are identied and the file is present at: GWAS_results/1_Final result.csv")),
        footer=tagList(actionButton("ok_gwas","Proceed to download the files"))
      ))

      observeEvent(input$ok_gwas, {
        removeModal()
      })

      output$man <- renderImage({
        manhattan_path <- file.path(dir, "GWAS_results", "1_Manhattan plot.jpeg")

        if (file.exists(manhattan_path)) {
          list(src = manhattan_path,
               contentType = 'image/jpeg',
               width = "100%",
               height = 350,
               alt = "Manhattan Plot")
        } else {
          return(NULL)  # Avoid errors if the file doesn't exist
        }
      }, deleteFile = FALSE)

      output$QQ <- renderImage({
        qqplot_path <- file.path(dir, "GWAS_results", "1_qq plot.jpeg")

        if (file.exists(qqplot_path)) {
          list(src = qqplot_path,
               contentType = 'image/jpeg',
               width = 400,
               height = 300,
               alt = "QQ Plot")
        } else {
          return(NULL)
        }
      }, deleteFile = FALSE)


      #we are chaning the above to all rs ids even if they come once
      res <- read.csv(file = file.path(dir,"GWAS_results","1_Final result.csv"))
      out <- res[!duplicated(res$RS.),]

      for (nr in c(1:nrow(out))) {
        out$Method[nr] <- toString(res$Method[res$RS.==out$RS.[nr]])
        out$QTN.effect[nr] <- mean(res$QTN.effect[res$RS.==out$RS.[nr]])
        out$r2....[nr] <- mean(res$r2....[res$RS.==out$RS.[nr]])
      }

      candgene_ind <- out[,c(5,6)]
      write_csv(candgene_ind, file.path(dir, "GWAS_results/GWAS_pos.csv"))
      write.csv(out[,c(2:8)],file.path(dir,"GWAS_results/Identified_MTAs.csv"),row.names = F)       ##-------- download

      output$result_table = renderTable({
        return(out)
      })

      output$downloadData <- downloadHandler(                                     ### Download GWAS results in Zip
        filename <- function() {
          paste("GAS_output", "zip", sep=".")
        },
        content <- function(file) {
          files2zip <- dir('GWAS_results', full.names = TRUE)
          zip(zipfile = 'testZip', files = files2zip)
          file.copy("testZip.zip", file)
        },
        contentType = "application/zip"
      )

    }else{
      output$analysis_complete = renderText({
        print("....GWAS analysis started....")
      })
    }
  })

  #Et-GWAS tabs----

  
  phe = reactive({
    phe <- read.csv(file = input$etgwas_pheno$datapath,header = T) #file name  =  pheno.csv
    
    pnam <- c("X.Phenotype.","trait") #Change trait
    colnames(phe) <- pnam
    return(phe)
  })
  
  plotData <- reactiveVal(NULL)
  processingResult <- reactiveVal(NULL)
  
  #dir.create(file.path(dir,"EtGWAS_results"))
  
  observeEvent(input$run_etgwas, {
    processingResult(NULL)
    
    showModal( modalDialog(
      h4(paste0("Et-GWAS for ",nrow(phe())," genotypes")),
      footer=tagList(h3("running..."))
    ))
    
    plotData(extract_irisID(trait = input$Trait,
                            infile = input$etgwas_pheno$datapath,
                            perc = as.numeric(input$Bulk_size),dir = dir,ip_dir = ip_dir))
    
    processingResult()
    
    et_out <- paste(input$Trait,"_intermediate_result",
                    as.numeric(input$Bulk_size),".csv",sep = "")
    
    et_fout <- paste(input$Trait,"_Final_result",
                     as.numeric(input$Bulk_size),".csv",sep = "")
    
    
    if(file.exists(file.path(dir,"EtGWAS_results",et_out))){
      removeModal()
      
      showModal( modalDialog(
        h4(paste0("MTAs are identied and the file is present at: EtGWAS_results/",
                  et_fout)),
        footer=tagList(actionButton("ok_etgwas","Proceed to download the files"))
      ))
      
      observeEvent(input$ok_etgwas, {
        removeModal()
      })
    }
  })
  
  output$resultText <- renderPrint({
    processingResult()
  })
  
  output$plot1 <- renderPlot({
    if (!is.null(plotData()))
      plotData()$plot5
  })
  output$plot2 <- renderPlot({
    if (!is.null(plotData()))
      plotData()$plot1
  })
  
  output$plot3 <- renderImage({
    outfile <- tempfile(fileext = '.jpeg')
    list(src = file.path(dir,"EtGWAS_results",
                         paste0(input$Trait,"_",
                                input$Bulk_size,"Manhattan.jpg")),
         contentType = 'image/jpeg',
         width = 400,
         height = 300,
         alt = "This is alternate text")
  }, deleteFile = F)
  
  output$downloadEt_Data <- downloadHandler(         ##### ### Download Et-GWAS results in Zip
    filename <- function() {
      paste("Et-output", "zip", sep=".")
    },
    content <- function(file) {
      files2zip <- dir('EtGWAS_results', full.names = TRUE)
      zip(zipfile = 'testZip', files = files2zip)
      file.copy("testZip", file)
    },
    contentType = "application/zip"
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Hap_pheno tabs----
  observeEvent(input$runhap, {

    # Set gff_path based on user selection
    if (input$choose_geno == "option1") {
      gff_path <- paste0(dir, "/Oryza_sativa.IRGSP-1.0.60.gff3")  # Default GFF
    } else if (input$choose_geno == "option2" && !is.null(input$gff_file)) {
      gff_path <- input$gff_file$datapath  # User-provided GFF
    } else {
      stop("Error: No GFF file provided!")
    }
    
    hap_path <- input$hapfile$datapath
    pos_path <- input$posfile$datapath
    ld_value <- input$LD
    select_criteria <- input$select_cri
    
    if (input$choose_geno == "option1") {
      dir <- getwd()
      gff_path <- paste0(dir, "/Oryza_sativa.IRGSP-1.0.60.gff3")  # Default GFF
    } else if (input$choose_geno == "option2" && !is.null(input$gff_file)) {
      gff_path <- input$gff_file$datapath  # User-provided GFF
    } else {
      stop("Error: No GFF file provided!")
    }
    
    # df <- candidate_gene(input$posfile$datapath, input$LD, gff_path)
    df <- candidate_gene(pos_path, ld_value, gff_path)

    # # df <- candidate_gene(pos_path, ld_value, input$gff_file$datapath)
    # df <- candidate_gene(pos_path, ld_value, gffpath)

    df <- as.data.frame(df)                                             ###############
    print(colnames(df))
    if (!"gene_id" %in% colnames(df)) {
      stop("Error: Column 'gene_id' not found in df.")
    }
    # df2 <- df[["gene_id"]]
    # df2 <- as.data.frame(df2)
    df2 <- df %>%
      dplyr::select(gene_id) %>%
      dplyr::distinct()
    write.csv(df2, file.path(dir, "locus.csv"), row.names = FALSE)

    showModal( modalDialog(
      h4(paste0("Haplo-Pheno analysis for ",nrow(df2)," Candidate genes")),
      footer=tagList(h3("running..."))
    ))

    dir1 <- file.path(getwd(), "Haplopheno")
    if (!dir.exists(dir1)) {
      dir.create(dir1, recursive = TRUE)
    }
    # dir <- getwd()
    if(file.exists(file.path(dir,"locus.csv"))){

      pheno <- read.csv(input$hapfile$datapath, header = TRUE)
      # pheno <- read.csv(hap_path, header = TRUE)
      all_pheno <- pheno[, 1:2]
      trait_dir <- colnames(all_pheno)[2] 
      # if (input$choose_geno == "option1") {
      #   all_pheno <- pheno
      # } else if (input$choose_geno == "option2") {
      #   if (ncol(pheno) < 2) {
      #     stop("Error: External phenotype file must have at least 2 columns.")         ################
      #   }
      #   all_pheno <- pheno[, 1:2]  # This overwrites the file path!
      # }

      # reactive object supplied with a forigen function
      gene <- hap_phe2(gene_infile = file.path(dir,"locus.csv"),
                       all_pheno = all_pheno,
                       select_cri = select_criteria,
                       dir1 = dir1)

      # func_piechart(gene_infile = file.path(dir,"locus.csv"),
      #               all_pheno = all_pheno,dir1 = dir1)
      
      observeEvent(c(input$season_pie, input$loc_id), {
        # req(input$season_pie, input$loc_id)
        req(input$loc_id)
        
        func_piechart(
          gene_infile = file.path(dir, "locus.csv"),
          all_pheno = all_pheno,
          dir1 = dir1
        )
      })
      
      
      removeModal()
    }

    shap <- reactive({
      haplotype_file <- file.path(dir1, "superior_haplotypes.csv")

      if (!file.exists(haplotype_file)) {
        return(data.frame(Message = "Error: Haplotype file not found"))
      }

      rfile <- read.csv(haplotype_file)
      rfile <- rfile[(!is.na(rfile$SH) & rfile$SH != "no"),]

      if (nrow(rfile) == 0) {
        return(data.frame(Message = "No superior haplotypes found"))
      }
      return(rfile)
    })

    showModal(modalDialog(
      h4(paste0("Haplo-Pheno analysis for is complete")),
      footer=tagList(actionButton("ok","Proceed to download the files"))
    ))

    observeEvent(input$ok, {
      removeModal()
    })

    output$Genes_table <- renderTable({
      req(file.exists("genes.csv"))
      head(read.csv("genes.csv"), 10)
    })

    output$table_down2 <- downloadHandler(
      filename = function() { "genes.csv" },
      content = function(file) {
        file.copy("genes.csv", file)
      }
    )

    output$candtable = renderTable({
      head(df)
    })

    output$canddw = downloadHandler(
      filename =  function(){
        "candidate_genes.csv"
      },
      content = function(file){
        write.csv(df,file)
      }
    )

    output$haptab = renderTable({
      head(gene)
    })

    output$hapdw = downloadHandler(
      file = function(){
        "superior_haplotype.csv"
      },
      content = function(file){
        write.csv(gene,file)
      }
    )

    warning_text <- reactive({
      if(nrow(shap()) == 0) {
        return(paste("<span style=\"color:red\">There are no superior haplotypes found</span>"))
      }
    })

    output$text1 <- renderText(warning_text())

    # updateSelectInput(session,"loc_id",choices = shap()[,1])
    updateSelectInput(session,"loc_id",choices = list.dirs(dir1, recursive = FALSE, full.names = FALSE))

    values <- reactiveValues()
    
    # output$piechart <- renderImage({
    #   # req(input$season_pie, input$loc_id)
    #   req(input$loc_id)
    #   
    #   # img_path <- file.path(dir1, input$season_pie, input$loc_id, "hap_subset_diversity.png")
    #   img_path <- file.path(dir1, input$loc_id, "hap_subset_diversity.png")
    #   message("Looking for image at: ", img_path)
    #   
    #   if (file.exists(img_path)) {
    #     list(
    #       src = img_path,
    #       contentType = 'image/png',
    #       width = 400,
    #       height = 300,
    #       alt = "Haplotype frequency pie chart"
    #     )
    #   } else {
    #     message("Image file not found at ", img_path)
    #     return(NULL)
    #   }
    # }, deleteFile = FALSE)
    
    output$piechart <- renderImage({
      req(input$loc_id)
      # filename <- file.path(dir1, input$loc_id, input$season_pie, "hap_subset_diversity.png")
      filename <- file.path(dir1, input$loc_id, trait_dir, "hap_subset_diversity.png")
      print(paste("Trying to load:", filename))  # debug line
      list(src = filename, contentType = 'image/png', width = 500, height = 500)
    }, deleteFile = FALSE)
    
    
    

    output$donortab <- renderTable({
      # req(input$loc_id, input$season_pie)
      req(input$loc_id)
      # donor_file <- file.path(dir1, input$loc_id, input$season_pie, paste0(input$loc_id, "_donors.csv"))
      donor_file <- file.path(dir1, input$loc_id, trait_dir, paste0(input$loc_id, "_donors.csv"))
      if (!file.exists(donor_file)) {
        return(data.frame(Message = "No donor file found for this locus."))
      }
      
      read.csv(donor_file) %>% head(5)
    })


    output$downloadhapData <- downloadHandler(
      filename <- function() {
        paste("output", "zip", sep=".")
      },
      content <- function(file) {
        dir.create("outhap")
        system("cp -r LOC* outhap")
        files2zip <- dir('outhap', full.names = TRUE)
        zip(zipfile = 'testZip', files = files2zip)
        file.copy("testZip.zip", file)
      },
      contentType = "application/zip"
    )
  })
}
