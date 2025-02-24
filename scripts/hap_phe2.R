hap_phe2 <- function(gene_infile, pheno_file, select_cri, dir1){
  dir1 <- file.path(getwd(), "Haplopheno")
  if (!dir.exists(dir1)) {
    dir.create(dir1, recursive = TRUE)
  }
  genes1<- read.csv(gene_infile, header = T)
  all_pheno <- read.csv(pheno_file)
  names(all_pheno)[1] <- "X.Phenotype."
  
  colnames(genes1) <- "gene_id"
  genes <- unique(genes1$gene_id)
  
  sup.hap <- matrix(nrow = 0,ncol = 3)
  sup.hap <- data.frame(sup.hap)
  colnames(sup.hap) <- c("LOC_ID","SH","Other.haps")
  
  #---------------------------------------------------------------------------
  duncan.multiple <-function (y, trt, DFerror, MSerror, alpha=0.05, group=TRUE,main = NULL,console=FALSE)
  {
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    if(is.null(main))main<-paste(name.y,"~", name.t)
    clase<-c("aov","lm")
    if("aov"%in%class(y) | "lm"%in%class(y)){
      if(is.null(main))main<-y$call
      A<-y$model
      DFerror<-df.residual(y)
      MSerror<-deviance(y)/DFerror
      y<-A[,1]
      ipch<-pmatch(trt,names(A))
      nipch<- length(ipch)
      for(i in 1:nipch){
        if (is.na(ipch[i]))
          return(if(console)cat("Name: ", trt, "\n", names(A)[-1], "\n"))
      }
      name.t<- names(A)[ipch][1]
      trt <- A[, ipch]
      if (nipch > 1){
        trt <- A[, ipch[1]]
        for(i in 2:nipch){
          name.t <- paste(name.t,names(A)[ipch][i],sep=":")
          trt <- paste(trt,A[,ipch[i]],sep=":")
        }}
      name.y <- names(A)[1]
    }
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    Mean<-mean(junto[,1])
    CV<-sqrt(MSerror)*100/Mean
    medians<-tapply.stat(junto[,1],junto[,2],stat="median")
    for(i in c(1,5,2:4)) {
      x <- tapply.stat(junto[,1],junto[,2],function(x)quantile(x)[i])
      medians<-cbind(medians,x[,2])
    }
    medians<-medians[,3:7]
    names(medians)<-c("Min","Max","Q25","Q50","Q75")	
    means <- tapply.stat(junto[,1],junto[,2],stat="mean") # change
    sds <-   tapply.stat(junto[,1],junto[,2],stat="sd") #change
    nn <-   tapply.stat(junto[,1],junto[,2],stat="length") # change
    means<-data.frame(means,std=sds[,2],r=nn[,2],medians)
    names(means)[1:2]<-c(name.t,name.y)
    ntr<-nrow(means)
    Tprob<-NULL
    k<-0
    for(i in 2:ntr){
      k<-k+1
      x <- suppressWarnings(warning(qtukey((1-alpha)^(i-1), i, DFerror)))
      Tprob[k]<-x
    }
    Tprob[Tprob == "NaN"] <- 0
    if(k<(ntr-1)){
      for(i in k:(ntr-1)){
        f <- Vectorize(function(x)ptukey(x,i+1,DFerror)-(1-alpha)^i)
        Tprob[i]<-uniroot(f, c(0,100))$root
      }
    }
    Tprob<-as.numeric(Tprob)
    nr <- unique(nn[,2])
    # Critical Value of Studentized Range
    if(console){
      cat("\nStudy:", main)
      cat("\n\nDuncan's new multiple range test\nfor",name.y,"\n")
      cat("\nMean Square Error: ",MSerror,"\n\n")
      cat(paste(name.t,",",sep="")," means\n\n")
      print(data.frame(row.names = means[,1], means[,2:6]))
    }
    if(length(nr) == 1 ) sdtdif <- sqrt(MSerror/nr)
    else {
      nr1 <-  1/mean(1/nn[,2])
      sdtdif <- sqrt(MSerror/nr1)
    }
    DUNCAN <- Tprob * sdtdif
    names(DUNCAN)<-2:ntr
    duncan<-data.frame(Table=Tprob,CriticalRange=DUNCAN)
    if ( group & length(nr) == 1 & console){
      cat("\nAlpha:",alpha,"; DF Error:",DFerror,"\n")
      cat("\nCritical Range\n")
      print(DUNCAN)
    }
    if ( group & length(nr) != 1 & console) cat("\nGroups according to probability of means differences and alpha level(",alpha,")\n")
    if ( length(nr) != 1) duncan<-NULL    
    Omeans<-order(means[,2],decreasing = TRUE) #correccion 2019, 1 abril.
    Ordindex<-order(Omeans)
    comb <-utils::combn(ntr,2)
    nn<-ncol(comb)
    dif<-rep(0,nn)
    DIF<-dif
    LCL<-dif
    UCL<-dif
    pvalue<-dif
    odif<-dif
    sig<-NULL
    for (k in 1:nn) {
      i<-comb[1,k]
      j<-comb[2,k]
      dif[k]<-means[i,2]-means[j,2]
      DIF[k]<-abs(dif[k])
      nx<-abs(i-j)+1
      odif[k] <- abs(Ordindex[i]- Ordindex[j])+1
      pvalue[k]<- round(1-ptukey(DIF[k]/sdtdif,odif[k],DFerror)^(1/(odif[k]-1)),4)
      LCL[k] <- dif[k] - DUNCAN[odif[k]-1]
      UCL[k] <- dif[k] + DUNCAN[odif[k]-1]
      sig[k]<-" "
      if (pvalue[k] <= 0.001) sig[k]<-"***"
      else  if (pvalue[k] <= 0.01) sig[k]<-"**"
      else  if (pvalue[k] <= 0.05) sig[k]<-"*"
      else  if (pvalue[k] <= 0.1) sig[k]<-"."
    }
    if(!group){  
      tr.i <- means[comb[1, ],1]
      tr.j <- means[comb[2, ],1]
      comparison<-data.frame("difference" = dif, pvalue=pvalue,"signif."=sig,LCL,UCL)
      rownames(comparison)<-paste(tr.i,tr.j,sep=" - ")
      if(console){cat("\nComparison between treatments means\n\n")
        print(comparison)}
      groups=NULL
    }
    if (group) {
      comparison=NULL
      # The probability matrix
      Q<-matrix(1,ncol=ntr,nrow=ntr)
      p<-pvalue
      k<-0
      for(i in 1:(ntr-1)){
        for(j in (i+1):ntr){
          k<-k+1
          Q[i,j]<-p[k]
          Q[j,i]<-p[k]
        }
      }
      groups <- orderPvalue(means[, 1], means[, 2],alpha, Q,console)
      names(groups)[1]<-name.y
      if(console) {
        cat("\nMeans with the same letter are not significantly different.\n\n")
        print(groups)
      }      
    }
    parameters<-data.frame(test="Duncan",name.t=name.t,ntr = ntr,alpha=alpha)
    statistics<-data.frame(MSerror=MSerror,Df=DFerror,Mean=Mean,CV=CV)
    rownames(parameters)<-" "
    rownames(statistics)<-" "
    rownames(means)<-means[,1]
    means<-means[,-1]
    output<-list(statistics=statistics,parameters=parameters, duncan=duncan,
                 means=means,comparison=comparison,groups=groups)
    class(output)<-"group"
    invisible(output)
  }
  ###--------------------- Superior Haplotype code -----------------------------
  for (gun in genes) {
    print(paste0(gun,": running"))
    dir.create(file.path(dir1,gun))
    colm <- colnames(all_pheno)
    colm<- colm[-1]
    file_path <- paste0("snpsift/",paste0(gun,".csv"))
    
    file.copy(from = file_path,to = paste0(dir1,"/",gun,"/",gun,".csv"))
    
    if(!file.exists(paste0(dir1,"/",gun,"/",gun,".csv"))){
      print(paste0("Candidate gene ",gun," doesn't exist. Skipping.."))
      next
    }
    
    gc<-c(gun)
    
    for (jun in colm) {
      dir.create(file.path(dir1,gun,jun))
      P <-all_pheno %>% dplyr::select(starts_with(jun))
      row.names(P) <- all_pheno$X.Phenotype.
      write.csv(P ,file.path(dir1,gun,"pheno.csv"),row.names = TRUE)
      
      haplo_P <- file.path(dir1,gun,jun)
      z<- haplo_P
      # setwd(paste0(z))
      
      filenames <- paste0(gun,".csv")
      filenames<- gsub(".csv","",filenames)
      
      y<- filenames
      a <- read.csv(paste0(dir1,"/",gun,"/",y,".csv"), header = F, stringsAsFactors = FALSE)
      snpfile<- y
      if (ncol(a) <=1) {
        #dir.create(paste0(snpfile,""))
        write.csv(a,file.path(dir1,gun,jun,"NULL.csv"))
        gcc<- "all"
        
      } 
      
      else if (ncol(a) <= 2) {
        
        gy<-t(a[1,])
        
        colnames(a)<-gy
        a<- data.frame(a[-c(1),])
        
        a[is.na(a)] <- "-"
        
        b <- a %>% filter(if_all(everything(), ~ !grepl("/",.)))
        b <- b %>% filter(if_all(everything(), ~ !grepl("N",.)))
        
        colnames(b)[1] <- c("ASSAY.ID")
        
        c<- gsub('.', '-', b$ASSAY.ID,fixed = T)
        
        c<- cbind(c,b[,-c(1:1)])
        
        c<-data.frame(c)
        
        colnames(c)<-colnames(b)
        
        yu<- c
        
        colnames(yu)<- c("A1","A2")
        
        as<- yu[!(str_count(yu$A2) < 1),]
        
        #dir.create(paste0(snpfile,""))
        
        sink(file.path(dir1,gun,jun,"Sequence.fas"))
        
        for (j in 1:nrow(as)){
          name = paste0(">",as[j,1])
          sequence = paste0(as[j,2])
          cat(name,sep="\n")
          cat(sequence,sep="\n")
        }
        sink()
        
        x<- read.fas(file.path(dir1,gun,jun,"Sequence.fas"))
        
        xi <- data.frame(x@sequence)
        
        haploy <- data.frame(table(xi$X1))
        
        colnames(haploy) <- c(gsub("X","",colnames(c)[-1]),"Counts")
        cnam <- gsub("X","",colnames(c)[-1])
        
        haploy$Freq <- haploy$Counts/nrow(as)
        
        popul <- paste("pop", seq(1:nrow(haploy)), sep = "")
        
        groups <- paste0("H",1:nrow(haploy))
        
        row.names(haploy) <- groups
        haploy$groups <- groups
        
        variety<- data.frame(as)
        
        var_pheno <- merge(variety,haploy, by.x = "A2",by.y = cnam)
        
        var_pheno <- data.frame(var_pheno [,c(1,2,5)])
        
        colnames(var_pheno) <- c("A2","A1","A3")
        
        var_pheno <-var_pheno %>% dplyr::select(order(colnames(var_pheno)))
        
        cnam2 <- colnames(c)[-1]
        
        colnames(var_pheno) <- c("ASSAY.ID",cnam2,"group")
        
        genotype <- data.frame(haploy[,1])
        
        row.names(genotype) <- groups
        
        colnames(genotype) <- colnames(c)[2]
        
        genotype[,1] <-as.character(genotype[,1]) 
        
        #genotype <- rbind(colnames(genotype), data.frame(genotype))
        
        #genotype <- cbind(row.names(genotype), data.frame(genotype))
        
        #genotype$`row.names(genotype)`[1]<-""
        
        write.csv(haploy[,-c(4)],file.path(dir1,gun,jun,"Haplotype.csv"))
        
        write.csv(var_pheno,file.path(dir1,gun,jun,"Variety.csv"))
        
        write.table(genotype,file.path(dir1,gun,jun,"Flapjack.genotype"),sep="\t", quote = F,row.names = F,col.names = F)
        
        hig<-str_remove(z,"/genes")
        
        pheno<- read.csv(paste0(dir1,"/",gun,"/pheno.csv"),header = T)
        
        pheno <- subset(pheno, !is.na(pheno[,2]))
        
        colnames(pheno)[1]<-"Designation"
        
        pheno$Designation <- gsub('IRIS ', 'IRIS-', pheno$Designation,fixed = T)
        
        haplopheno <- merge(var_pheno, pheno, by.x = 'ASSAY.ID', by.y = 'Designation')
        
        # haplopheno$group<- as.factor(haplopheno$group)
        
        haplopheno$group<- as.character(haplopheno$group)
        
        haplopheno<-haplopheno[order (haplopheno$group),]
        
        hp <- plyr::count(haplopheno$group)
        
        hp <- hp[!(hp$freq <3),]
        
        if(length(hp[,1]) >=2)
        {
          sd<- merge(hp,haplopheno, by.x = 'x', by.y = 'group')
          
          sd<- sd  %>% dplyr::select(1,ncol(sd))
          
          colu<- colnames(sd)
          
          colnames(sd)<- c("Haplotype",colu[2])
          
          colu <- colnames(sd)
          
          sd<- subset(sd, !is.na(sd[,2]))
          
          sink(file.path(dir1,gun,jun,"anova.txt"))
          
          aooov<- aov(formula((paste0(colnames(sd)[2], "~", colnames(sd)[1]))), data = sd)
          
          print(summary(aooov))
          
          sink()
          
          if( sum(aooov$effects) == 0){
            sink(file.path(dir1,gun,jun, "no_anova.txt"))
            print(" no duncan test can be performed contrasts can be applied only to factors with 2 or more levels")
            sink()
            gcc<-"NA"
          }
          else {
            
            sink(file.path(dir1,gun,jun,"duncan.txt"))
            print(duncan<- duncan.multiple(aooov,trt = "Haplotype",console = TRUE, group = T))
            sink()
            
            sink(file.path(dir1,gun,jun,"duncan_pval.txt"))
            
            print(duncan.pval<- duncan.multiple(aooov,trt = "Haplotype",console = TRUE, group = F))
            
            sink()
            
            summary<- data.frame(duncan$means)
            
            #summary<- cbind(hp,summary)
            
            write.csv(summary,file.path(dir1,gun,jun,"summary.stat.csv"))
            
            labels <- data.frame(duncan$groups)
            
            labels$Haplotype<- row.names(labels)
            
            labels<- labels[order(labels$Haplotype),]
            
            hp <- subset(hp,hp$x%in% labels$Haplotype)
            
            specify<- cbind(hp,labels)
            
            mycolors = c(brewer.pal(name="BuGn", n = 9), brewer.pal(name="OrRd", n = 9),brewer.pal(name = "YlOrBr", n=9),brewer.pal(name = "Pastel1",n=9))
            
            jpeg(file.path(dir1,gun,jun,"boxplot.jpg"), height = 1200,width = 1400, res = 300)
            
            p<-ggplot(sd, aes(x=Haplotype, y=sd[,2])) +labs(y=paste0(colnames(sd)[2]))+ geom_boxplot(aes(fill = factor(Haplotype)), show.legend = F)+scale_color_manual(values = mycolors, aesthetics = "fill")+ theme(axis.text.x = element_text(angle = 90, hjust= 1.0))
            
            q<-p+geom_text(data = specify, aes(x,Inf, label= paste0("n=",freq,"\n",groups)), vjust = 1,size= 1.5)
            
            r<- q+theme_bw()+theme(text = element_text(size=10, face="bold"),axis.text.x = element_text(angle = 90, hjust= 1.0))
            
            s<- r + geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5,color = "red")
            
            t <- s+ggtitle(label = paste0(snpfile,""))
            
            print(t)
            
            dev.off()
            
            var_pop <- gsub('_', ' ', var_pheno[,1])
            
            variety_pop <- var_pheno
            
            variety_pop$names <- var_pop
            
            variety_pop<- merge(variety_pop, pheno, by.x = 'names', by.y = 'Designation')
            
            variety_pop<- variety_pop[order (variety_pop$group),]
            
            variety_pop <- data.frame(variety_pop [,2:5])
            
            colnames(variety_pop)[1] <- "Designation"
            
            write.csv(variety_pop,file.path(dir1,gun,jun,"variety_subset.csv"),row.names = F)
            
            count <- plyr::count(haplopheno$group)
            
            colnames(count)<- c("x", "count_pop")
            
            count$freq_pop <- count$count_pop/nrow(haplopheno)
            
            hap<- cbind(row.names(haploy), haploy)
            
            haplotype_pop<- merge(hap, count, by.x = 'row.names(haploy)', by.y = 'x' )
            
            haplotype_pop<- data.frame(haplotype_pop[,-c(5)])
            
            colnames(haplotype_pop)[1:2] <- c("Haplotype",gsub("X","",colnames(haplotype_pop)[2]))
            
            write.csv(haplotype_pop,file.path(dir1,gun,jun,"haplotype_subset.csv"),row.names = F)
            
            gcc<- "all"
            
          }
          
        }  
        
        else {
          sink(file.path(dir1,gun,jun,"no_anova.txt"))
          print(" no duncan test can be performed contrasts can be applied only to factors with 2 or more levels")
          sink()
          gcc<- "NA"
        }
        
        
      }
      else {
        gy<-t(a[1,])
        
        colnames(a)<-gy
        
        a<- data.frame(a[-c(1),])
        
        a[is.na(a)] <- "-"
        
        a[is.null(a)] <- "-"
        
        b <- a %>% filter(across(everything(a), ~ !grepl("/",.)))
        b <- b %>% filter(across(everything(b), ~ !grepl("N",.)))

        colnames(b)[1]<- c("ASSAY.ID")
        
        c<- gsub('.', '-', b$ASSAY.ID,fixed = T)
        
        c<- cbind(c,b[,-c(1:1)])
        
        c<-data.frame(c)
        
        colnames(c)<-colnames(b)
        
        d<- unite(c[,-c(1:1)], col = "seq", remove =  TRUE, sep = "", na.rm = TRUE)
        
        e <- data.frame(cbind(c$ASSAY.ID,d$seq))
        
        e<-e[!(str_count(e$X2) < 1),]
        
        #dir.create(paste0(snpfile,""))
        
        #sink(file.path((paste0(snpfile,"")), "Sequence.fas"))
        sink(file.path(dir1,gun,jun,"Sequence.fas"))
        
        for (j in 1:nrow(e)){
          name = paste0(">",e[j,1])
          sequence = paste0(e[j,2])
          cat(name,sep="\n")
          cat(sequence,sep="\n")
        }
        sink()
        
        f<- data.frame(table(e$X2))
        
        f$group <- paste0("H",seq(1:nrow(f)))
        
        colnames(f)[1] <- "X2"
        
        g<-join(e,f,by = 'X2')
        
        variety <- cbind(c,g[4])
        
        haplo <- data.frame(strsplit(as.character(f$X2),split = ""))
        
        haplo<- data.matrix(t(haplo))
        
        row.names(haplo) <- f$group
        
        colnames(haplo) <- gsub("X","",colnames(c)[-1])
        
        haplo <- cbind(haplo,f[2])
        
        colnames(haplo)[ncol(haplo)] <- "Counts"
        
        haplo$Freq <- haplo$Counts/nrow(e)
        
        var_pheno <- variety
        
        colnames(var_pheno)[1] <- c("Accession")
        
        genotype <- data.frame(rev((haplo))[-c(1:2)])
        
        genotype <- data.frame(rev((genotype)))
        
        #variety<- merge(variety,c, by.x = 'e...1.1.',by.y = 'ASSAY.ID')
        
        write.csv(haplo,file.path(dir1,gun,jun,"Haplotype.csv"))
        
        write.csv(variety,file.path(dir1,gun,jun,"Variety.csv"))
        
        write.table(genotype,file.path(dir1,gun,jun,"Flapjack.genotype"),sep="\t", quote = F,row.names = F,col.names = F)
        
        #hig<-str_remove(z,"/genes")
        
        hig<-str_remove(z,"/genes")
        
        pheno<- read.csv(paste0(dir1,"/",gun,"/","pheno.csv"),header = T)
        
        pheno <- subset(pheno, !is.na(pheno[,2]))
        
        colnames(pheno)[1]<-"Designation"
        
        pheno$Designation <- gsub('IRIS ', 'IRIS-', pheno$Designation,fixed = T)
        
        haplopheno <- merge(var_pheno, pheno, by.x = 'Accession', by.y = 'Designation')
        
        # haplopheno$group<- as.factor(haplopheno$group)
        
        haplopheno$group<- as.character(haplopheno$group)
        
        haplopheno<-haplopheno[order (haplopheno$group),]
        
        hp <- plyr::count(haplopheno$group)
        
        hp <- hp[!(hp$freq <3),]
        
        if(length(hp[,1])>=2)
        {
          sd<- merge(hp,haplopheno, by.x = 'x', by.y = 'group')
          
          sd<- sd  %>% dplyr::select(1,ncol(sd))
          
          colu<- colnames(sd)
          
          colnames(sd)<- c("Haplotype",colu[2])
          
          colu <- colnames(sd)
          
          sd<- subset(sd, !is.na(sd[,2]))
          
          sink(file.path(dir1,gun,jun,"anova.txt"))
          
          aooov<- lm(formula((paste0(colnames(sd)[2], "~", colnames(sd)[1]))), data = sd)
          
          print(summary(aooov))
          
          sink()
          
          if( sum(aooov$effects) == 0){
            sink(file.path(dir1,gun,jun, "no_anova.txt"))
            print(" no duncan test can be performed contrasts can be applied only to factors with 2 or more levels")
            sink()
            gcc<-"NA"
          }
          else {
            
            sink(file.path(dir1,gun,jun,"duncan.txt"))
            
            print(duncan<- duncan.multiple(aooov,trt = "Haplotype",console = FALSE, group = TRUE,main = TRUE))
            
            sink()
            
            sink(file.path(dir1,gun,jun,"duncan_pval.txt"))
            
            print(duncan.pval<- duncan.multiple(aooov,trt = "Haplotype",console = TRUE, group = F))
            
            sink()
            
            summary<- data.frame(duncan$means)
            
            #summary<- cbind(hp,summary)
            
            write.csv(summary,file.path(dir1,gun,jun,"summary.stat.csv"))
            
            labels <- data.frame(duncan$groups)
            
            labels$Haplotype<- row.names(labels)
            
            labels<- labels[order(labels$Haplotype),]
            
            hp <- subset(hp,hp$x%in% labels$Haplotype)
            
            specify<- cbind(hp,labels)
            
            mycolors = c(brewer.pal(name="BuGn", n = 9), brewer.pal(name="OrRd", n = 9),brewer.pal(name = "YlOrBr", n=9),brewer.pal(name = "Pastel1",n=9),brewer.pal(name = "Set1",n=9),brewer.pal(name="BuGn", n = 9), brewer.pal(name="OrRd", n = 9),brewer.pal(name = "YlOrBr", n=9),brewer.pal(name = "Pastel1",n=9),brewer.pal(name = "Set1",n=9))
            
            jpeg(file.path(dir1,gun,jun,"boxplot.jpg"), height = 1200,width = 1400, res = 300)
            
            p<-ggplot(sd, aes(x=Haplotype, y=sd[,2])) +labs(y=paste0(colnames(sd)[2]))+ geom_boxplot(aes(fill = factor(Haplotype)), show.legend = F)+scale_color_manual(values = mycolors, aesthetics = "fill")+ theme(axis.text.x = element_text(angle = 90, hjust= 1.0))
            
            q<-p+geom_text(data = specify, aes(x,Inf, label= paste0("n=",freq,"\n",groups)), vjust = 1,size= 1.5)
            
            r<- q+theme_bw()+theme(text = element_text(size=10, face="bold"),axis.text.x = element_text(angle = 90, hjust= 1.0))
            
            s<- r + geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5,color = "red")
            
            t<- s+ggtitle(label = paste0(snpfile,""))
            
            print(t)
            
            dev.off()
            
            var_pop <- gsub('_', ' ', variety[,1])
            
            variety_pop <- variety
            
            variety_pop$names <- var_pop
            
            variety_pop<- merge(variety_pop, pheno, by.x = 'names', by.y = 'Designation')
            
            variety_pop<- variety_pop[order (variety_pop$group),]
            
            variety_pop <- data.frame(variety_pop [,-c(1:1)])
            
            colnames(variety_pop)[1] <- "Designation"
            
            write.csv(variety_pop,file.path(dir1,gun,jun,"variety_subset.csv"),row.names = F)
            
            count <- plyr::count(haplopheno$group)
            
            colnames(count)<- c("x", "count_pop")
            
            count$freq_pop <- count$count_pop/nrow(haplopheno)
            
            hap<- cbind(row.names(haplo), haplo)
            
            haplotype_pop<- merge(hap, count, by.x = 'row.names(haplo)', by.y = 'x' )
            
            colnames(haplotype_pop)[1] <- c("Haplotype")
            
            write.csv(haplotype_pop,file.path(dir1,gun,jun,"haplotype_subset.csv"),row.names = F)
            gcc<- "all"
          }
        }
        
        else {
          sink(file.path(dir1,gun,jun,"no_anova.txt"))
          print(" no duncan test can be performed contrasts can be applied only to factors with 2 or more levels")
          sink()
          gcc<- "NA"
        }
      }
      if( gcc=="NA"){
        #gcc<- "NA"
        gc<- c(gc,gcc,"NA")
        rm(gcc)
        #rm(aooov)
      }
      else  if (nrow(table(duncan$groups$groups))>=2){
        bea <- duncan$groups #object name changed
        come <- data.frame(table(duncan$groups$groups))
        
        #------------------------- Change
        if (select_cri == "High"){
          funi <- come[1:1,]
          funi = data.frame(funi)
        }else{
          funi <- come[nrow(come),]
          funi = data.frame(funi)
        }
        #--------------------------     
        happi <- funi$Var1
        
        happi <- subset(bea,bea$groups %in% happi)
        happi <- rownames(happi)
        hc <- happi
        
        mean_val <- subset(summary,rownames(summary) %in% happi)  # change
        other_mean_val <- subset(summary, !(rownames(summary) %in% happi))
        
        shap_mean <- str_c(row.names(mean_val),
                       paste("(", round(mean_val[,1],2),"±",round(mean_val[,2],2), ")", sep = ""),
                       collapse = ",") # change
        
        otherhap_mean <- str_c(row.names(other_mean_val),
                             paste("(", round(other_mean_val[,1],2),"±",round(other_mean_val[,2],2), ")", sep = ""),
                             collapse = ",") # change
        
        gc<-c(gc,shap_mean,otherhap_mean)
  
        donors <- data.frame()
        for (hun in 1:length(hc)) {
          tuguu <- hc[hun]
          fhi <- subset(variety_pop, variety_pop$group %in% tuguu)
          fhi <- fhi[,c(1,ncol(fhi))]
          fhi$SH <- tuguu
          pun1 <- subset(summary,rownames(summary) %in% tuguu)
          fhi$meanSH <- round(pun1[,1],digits = 2)
        }
        donors <- rbind(donors,fhi)
        write.csv(donors, file.path(dir1,gun,jun,paste0(gun,"_donors.csv")), 
                  row.names = F)
      }
      else {
        hf<- "no"
        allhap_mean <- str_c(row.names(summary),
                           paste("(", round(summary[,1],2),"±",round(summary[,2],2), ")", sep = ""),
                           collapse = ",") # change
        gc<- c(gc,hf, allhap_mean)
        #rm(aooov)
        #rm(duncan)
        #rm(duncan.pval)
      }
      rm(aooov)
      rm(duncan)
      rm(duncan.pval)
      # setwd(file.path(dir,gun))
      
    }
    #setwd(file.path(dir,gun))
  
    gc <- t(data.frame(gc))
    colnames(gc) <- c("LOC_ID","SH","Other.haps")
    sup.hap<-rbind(sup.hap,gc)
    rownames(sup.hap) <- NULL
    # setwd(dir)
    
  }
  colnames(sup.hap) <- c("LOC_ID","SH","Other.haps")
  
  write.csv(sup.hap,file = file.path(dir1,"superior_haplotypes.csv"),
            row.names = F)
  
  return(sup.hap)
}

