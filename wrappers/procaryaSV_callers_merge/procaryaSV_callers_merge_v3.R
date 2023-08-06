
library(data.table)
library(stringr)
library(seqinr)
library(ggVennDiagram)
library(ggplot2)
library(RColorBrewer)


run_all <- function(args){
  # arguments
  output <- args[1]
  venn_png <- args[2]
  barchart_png <- args[3]
  fastaFile <- args[4] 
  min_sv_length <- as.integer(args[5]) #default 1 
  max_sv_length <- as.integer(args[6]) #default refLength/5
  distanceThreshold <- as.integer(args[7])
  vcf_files <- args[8:length(args)]
  
  ## DEFAULT VALUES
  numCallers <- length(vcf_files)
  if(is.na(distanceThreshold)){distanceThreshold <- 100} #default
  if(is.na(min_sv_length)){min_sv_length <- 50} #default
  # max_sv_length depends on chromosome length
  
  ## read FASTA
  referenceLengthList <- seqinr::getLength(seqinr::read.fasta(fastaFile, seqonly=TRUE))
  numberOfFastaContigs <- length(referenceLengthList)
  fastaContiglist <- seqinr::getName(seqinr::read.fasta(fastaFile, seqonly=FALSE))
  
  ## detect files by caller
  delly <- vcf_files[which(grepl( "delly2", vcf_files, fixed = TRUE))]
  lumpy <- vcf_files[which(grepl( "lumpy", vcf_files, fixed = TRUE))]
  cnvnator <- vcf_files[which(grepl( "cnvnator", vcf_files, fixed = TRUE))]
  pindel <- vcf_files[which(grepl( "pindel", vcf_files, fixed = TRUE))]
  cnproscan <- vcf_files[which(grepl( "cnproscan", vcf_files, fixed = TRUE))]
  
  ##############################################################################
  ## LOAD VCF FILES
  
  #DELLY2
  if(!identical(delly, character(0))){
    dellyDF <- fread(file=delly, sep='\t', header = TRUE, skip = '#CHROM')
    if(nrow(dellyDF)){
      dellyDF <- dellyDF[!grepl("SVTYPE=BND", dellyDF$INFO),]
      dellyDF$END <- sapply(dellyDF$INFO, function(x)   str_extract(str_extract(x,"END=*?([0-9]+)"),"[0-9]+"))
      dellyDF$SVTYPE <- dellyDF$ALT
      dellyDF$CALLER <- "delly2"
      dellyDF$VCF_INFO <- paste0("QUAL:",dellyDF$QUAL,"\t","FILTER:",dellyDF$FILTER,"\t","INFO:",dellyDF$INFO)
      dellyDF <- dellyDF[,c("#CHROM","POS","END","SVTYPE","CALLER","VCF_INFO")]
      dellyDF$END <- as.numeric(dellyDF$END)
      dellyDF$POS <- as.numeric(dellyDF$POS)
      dellyDF$LEN <- dellyDF$END - dellyDF$POS
      dellyDF$SVTYPE <- gsub("<", "", dellyDF$SVTYPE) #remove < from SVTYPE
      dellyDF$SVTYPE <- gsub(">", "", dellyDF$SVTYPE) #remove < from SVTYPE
    }else{delly <- character(0)}
  }else{
    dellyDF <- data.table("#CHROM"=character(),"POS"=numeric(),"END"=numeric(),"SVTYPE"=character(),"CALLER"=character(),"VCF_INFO"=character(),"LEN"=numeric() )
  }
  
  # LUMPY
  if(!identical(lumpy, character(0))){
    lumpyDF <- fread(file=lumpy, sep='\t', header = TRUE, skip = '#CHROM')
    if(nrow(lumpyDF)){
      lumpyDF <- lumpyDF[!grepl("SVTYPE=BND", lumpyDF$INFO),]
      lumpyDF$SVTYPE <- lumpyDF$ALT
      lumpyDF$END <- sapply(lumpyDF$INFO, function(x)   str_extract(str_extract(x,"END=*?([0-9]+)"),"[0-9]+"))
      lumpyDF$CALLER <- "lumpy"
      lumpyDF$VCF_INFO <- paste0("INFO:",lumpyDF$INFO)
      lumpyDF <- lumpyDF[,c("#CHROM","POS","END","SVTYPE","CALLER","VCF_INFO")]
      lumpyDF$END <- as.numeric(lumpyDF$END)
      lumpyDF$POS <- as.numeric(lumpyDF$POS)
      lumpyDF$LEN <- lumpyDF$END - lumpyDF$POS
      lumpyDF$SVTYPE <- gsub("<", "", lumpyDF$SVTYPE) #remove < from SVTYPE
      lumpyDF$SVTYPE <- gsub(">", "", lumpyDF$SVTYPE) #remove < from SVTYPE
    }else{lumpy <- character(0)}
  }else{
    lumpyDF <- data.table("#CHROM"=character(),"POS"=numeric(),"END"=numeric(),"SVTYPE"=character(),"CALLER"=character(),"VCF_INFO"=character(),"LEN"=numeric() )
  }
  
  # cnvnator
  if(!identical(cnvnator, character(0))){
    cnvnatorDF <- fread(file=cnvnator, sep='\t', header = TRUE, skip = '#CHROM')
    if(nrow(cnvnatorDF)){
      cnvnatorDF <- cnvnatorDF[!grepl("BND", cnvnatorDF$INFO),]
      cnvnatorDF$SVTYPE <- cnvnatorDF$ALT
      cnvnatorDF$END <- sapply(cnvnatorDF$INFO, function(x)   str_extract(str_extract(x,"END=*?([0-9]+)"),"[0-9]+"))
      cnvnatorDF$CALLER <- "cnvnator"
      cnvnatorDF$VCF_INFO <- paste0("INFO:",cnvnatorDF$INFO)
      cnvnatorDF <- cnvnatorDF[,c("#CHROM","POS","END","SVTYPE","CALLER","VCF_INFO")]
      cnvnatorDF$END <- as.numeric(cnvnatorDF$END)
      cnvnatorDF$POS <- as.numeric(cnvnatorDF$POS)
      cnvnatorDF$LEN <- cnvnatorDF$END - cnvnatorDF$POS
      cnvnatorDF$SVTYPE <- gsub("<", "", cnvnatorDF$SVTYPE) #remove < from SVTYPE
      cnvnatorDF$SVTYPE <- gsub(">", "", cnvnatorDF$SVTYPE) #remove < from SVTYPE
    }else{cnvnator <- character(0)}
  }else{
    cnvnatorDF <- data.table("#CHROM"=character(),"POS"=numeric(),"END"=numeric(),"SVTYPE"=character(),"CALLER"=character(),"VCF_INFO"=character(),"LEN"=numeric() )
  }
  
  # pindel
  if(!identical(pindel, character(0))){
    pindelDF <- fread(file=pindel, sep='\t', header = TRUE, skip = '#CHROM')
    if(nrow(pindelDF)){
      pindelDF <- pindelDF[!grepl("BND", pindelDF$INFO),]
      pindelDF$END <- sapply(pindelDF$INFO, function(x)   str_extract(str_extract(x,"END=*?([0-9]+)"),"[0-9]+"))
      pindelDF$SVTYPE <- sapply(pindelDF$INFO, function(x)   gsub("SVTYPE=","",str_extract(x,"SVTYPE=[A-Z]+")))
      pindelDF$CALLER <- "pindel"
      pindelDF$VCF_INFO <- paste0("INFO:",pindelDF$INFO)
      pindelDF <- pindelDF[,c("#CHROM","POS","END","SVTYPE","CALLER","VCF_INFO")]
      pindelDF$END <- as.numeric(pindelDF$END)
      pindelDF$POS <- as.numeric(pindelDF$POS)
      pindelDF$LEN <- pindelDF$END - pindelDF$POS
      pindelDF$SVTYPE <- gsub("<", "", pindelDF$SVTYPE) #remove < from SVTYPE
      pindelDF$SVTYPE <- gsub(">", "", pindelDF$SVTYPE) #remove < from SVTYPE
    }else{pindel <- character(0)}
  }else{
    pindelDF <- data.table("#CHROM"=character(),"POS"=numeric(),"END"=numeric(),"SVTYPE"=character(),"CALLER"=character(),"VCF_INFO"=character(),"LEN"=numeric() )
  }
  
  # cnproscan
  if(!identical(cnproscan, character(0))){
    cnproscanDF <- fread(file=cnproscan, sep='\t', header = TRUE, skip = '#CHROM')
    if(nrow(cnproscanDF)){
      cnproscanDF <- cnproscanDF[!grepl("BND", cnproscanDF$INFO),]
      cnproscanDF$LEN <- sapply(cnproscanDF$INFO, function(x)   str_extract(str_extract(x,"SVLEN=*?([0-9]+)"),"[0-9]+"))
      cnproscanDF$END <- sapply(cnproscanDF$INFO, function(x)   str_extract(str_extract(x,"END=*?([0-9]+)"),"[0-9]+"))
      cnproscanDF$SVTYPE <- cnproscanDF$ALT
      cnproscanDF$CALLER <- "cnproscan"
      cnproscanDF$VCF_INFO <- paste0("INFO:",cnproscanDF$INFO)
      cnproscanDF <- cnproscanDF[,c("#CHROM","POS","END","SVTYPE","CALLER","VCF_INFO")]
      cnproscanDF$END <- as.numeric(cnproscanDF$END)
      cnproscanDF$POS <- as.numeric(cnproscanDF$POS)
      cnproscanDF$LEN <- cnproscanDF$END - cnproscanDF$POS
      cnproscanDF$SVTYPE <- gsub("<", "", cnproscanDF$SVTYPE) #remove < from SVTYPE
      cnproscanDF$SVTYPE <- gsub(">", "", cnproscanDF$SVTYPE) #remove < from SVTYPE
    }else{cnproscan <- character(0)}
  }else{
    cnproscanDF <- data.table("#CHROM"=character(),"POS"=numeric(),"END"=numeric(),"SVTYPE"=character(),"CALLER"=character(),"VCF_INFO"=character(),"LEN"=numeric() )
  }
  
  ############################################################################################################################################################
  ############################################################################################################################################################
  ## GO OVER CHROMOSOMES
  
  SHORT_RESULTS <- data.table()
  
  for (contigNum in 1:numberOfFastaContigs){
    
    refHeader <- fastaContiglist[contigNum]
    refLength <- referenceLengthList[contigNum]
    if(is.na(max_sv_length)){max_sv_length <- refLength/5} #default
    
    # DELLY2
    if(!identical(delly, character(0))){
      TMPdellyDF <- dellyDF[dellyDF$`#CHROM`==refHeader,]
      TMPdellyDF <- TMPdellyDF[TMPdellyDF$LEN<=max_sv_length,]  #
      TMPdellyDF <- TMPdellyDF[TMPdellyDF$LEN>=min_sv_length,]
      dellyVEC_DUP <- rep(0,refLength)
      dellyVEC_DEL <- rep(0,refLength)
      dellyVEC_INS <- rep(0,refLength)
      dellyVEC_INV <- rep(0,refLength)
      if(nrow(TMPdellyDF)){
        for (j in 1:nrow(TMPdellyDF)){
          if(TMPdellyDF$SVTYPE[j]=="DUP"){dellyVEC_DUP[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
          if(TMPdellyDF$SVTYPE[j]=="DEL"){dellyVEC_DEL[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
          if(TMPdellyDF$SVTYPE[j]=="INS"){dellyVEC_INS[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
          if(TMPdellyDF$SVTYPE[j]=="INV"){dellyVEC_INV[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
        }}
    }else{
      dellyVEC_DUP <- rep(0,refLength)
      dellyVEC_DEL <- rep(0,refLength)
      dellyVEC_INS <- rep(0,refLength)
      dellyVEC_INV <- rep(0,refLength)
      TMPdellyDF <- dellyDF
    }
    
    # LUMPY
    if(!identical(lumpy, character(0))){
      TMPlumpyDF <- lumpyDF[lumpyDF$`#CHROM`==refHeader,]
      TMPlumpyDF <- TMPlumpyDF[TMPlumpyDF$LEN<=max_sv_length,]  #
      TMPlumpyDF <- TMPlumpyDF[TMPlumpyDF$LEN>=min_sv_length,]
      lumpyVEC_DUP <- rep(0,refLength)
      lumpyVEC_DEL <- rep(0,refLength)
      if(nrow(TMPlumpyDF)){
        for (j in 1:nrow(TMPlumpyDF)){
          if(TMPlumpyDF$SVTYPE[j]=="DUP"){lumpyVEC_DUP[TMPlumpyDF$POS[j]:TMPlumpyDF$END[j]] <- 1}
          if(TMPlumpyDF$SVTYPE[j]=="DEL"){lumpyVEC_DEL[TMPlumpyDF$POS[j]:TMPlumpyDF$END[j]] <- 1}
        }}
    }else{
      lumpyVEC_DUP <- rep(0,refLength)
      lumpyVEC_DEL <- rep(0,refLength)
      TMPlumpyDF <- lumpyDF
    }
    
    # cnvnator
    if(!identical(cnvnator, character(0))){
      TMPcnvnatorDF <- cnvnatorDF[cnvnatorDF$`#CHROM`==refHeader,]
      TMPcnvnatorDF <- TMPcnvnatorDF[TMPcnvnatorDF$LEN<=max_sv_length,]  #
      TMPcnvnatorDF <- TMPcnvnatorDF[TMPcnvnatorDF$LEN>=min_sv_length,]
      cnvnatorVEC_DUP <- rep(0,refLength)
      cnvnatorVEC_DEL <- rep(0,refLength)
      if(nrow(TMPcnvnatorDF)){
        for (j in 1:nrow(TMPcnvnatorDF)){
          if(TMPcnvnatorDF$SVTYPE[j]=="DUP"){cnvnatorVEC_DUP[TMPcnvnatorDF$POS[j]:TMPcnvnatorDF$END[j]] <- 1}
          if(TMPcnvnatorDF$SVTYPE[j]=="DEL"){cnvnatorVEC_DEL[TMPcnvnatorDF$POS[j]:TMPcnvnatorDF$END[j]] <- 1}
        }}
    }else{
      cnvnatorVEC_DUP <- rep(0,refLength)
      cnvnatorVEC_DEL <- rep(0,refLength)
      TMPcnvnatorDF <- cnvnatorDF
    }
    
    # pindel
    if(!identical(pindel, character(0))){
      TMPpindelDF <- pindelDF[pindelDF$`#CHROM`==refHeader,]
      TMPpindelDF <- TMPpindelDF[TMPpindelDF$LEN<=max_sv_length,]  #
      TMPpindelDF <- TMPpindelDF[TMPpindelDF$LEN>=min_sv_length,]
      pindelVEC_DUP <- rep(0,refLength)
      pindelVEC_DEL <- rep(0,refLength)
      pindelVEC_INS <- rep(0,refLength)
      pindelVEC_INV <- rep(0,refLength)
      # ingoring RPL - replacement https://www.biostars.org/p/113765/
      if(nrow(TMPpindelDF)){
        for (j in 1:nrow(TMPpindelDF)){
          if(TMPpindelDF$SVTYPE[j]=="DUP"){pindelVEC_DUP[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
          if(TMPpindelDF$SVTYPE[j]=="DEL"){pindelVEC_DEL[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
          if(TMPpindelDF$SVTYPE[j]=="INS"){pindelVEC_INS[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
          if(TMPpindelDF$SVTYPE[j]=="INV"){pindelVEC_INV[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
        }}
    }else{
      pindelVEC_DUP <- rep(0,refLength)
      pindelVEC_DEL <- rep(0,refLength)
      pindelVEC_INS <- rep(0,refLength)
      pindelVEC_INV <- rep(0,refLength)
      TMPpindelDF <- pindelDF
    }
    
    # cnproscan
    if(!identical(cnproscan, character(0))){
      TMPcnproscanDF <- cnproscanDF[cnproscanDF$`#CHROM`==refHeader,]
      TMPcnproscanDF <- TMPcnproscanDF[TMPcnproscanDF$LEN<=max_sv_length,]  #
      TMPcnproscanDF <- TMPcnproscanDF[TMPcnproscanDF$LEN>=min_sv_length,]
      cnproscanVEC_DUP <- rep(0,refLength)
      cnproscanVEC_DEL <- rep(0,refLength)
      # ingoring RPL - replacement https://www.biostars.org/p/113765/
      if(nrow(TMPcnproscanDF)){
        for (j in 1:nrow(TMPcnproscanDF)){
          if(TMPcnproscanDF$SVTYPE[j]=="DUP"){cnproscanVEC_DUP[TMPcnproscanDF$POS[j]:TMPcnproscanDF$END[j]] <- 1}
          if(TMPcnproscanDF$SVTYPE[j]=="DEL"){cnproscanVEC_DEL[TMPcnproscanDF$POS[j]:TMPcnproscanDF$END[j]] <- 1}
        }}
    }else{
      cnproscanVEC_DUP <- rep(0,refLength)
      cnproscanVEC_DEL <- rep(0,refLength)
      TMPcnproscanDF <- cnproscanDF
    }
    
    
    ## CREATE VECTORS
    DELETIONS <- dellyVEC_DEL+lumpyVEC_DEL+cnvnatorVEC_DEL+pindelVEC_DEL+cnproscanVEC_DEL
    DUPLICATIONS <- dellyVEC_DUP+lumpyVEC_DUP+cnvnatorVEC_DUP+pindelVEC_DUP+cnproscanVEC_DUP
    INSERTIONS <- pindelVEC_INS+dellyVEC_INS
    INVERSIONS <- pindelVEC_INV+dellyVEC_INV
    
    ##############################################################################
    ## MERGE CLOSE SVs - #iterative for every level separately
    ## CLOSE GAPS IN BETWEEN TIGHT SVs DELETIONS
    if(as.integer(max(DELETIONS))>0){
      for(level in as.integer(max(DELETIONS)):1){
        nonzeroIDX <- sort(which(DELETIONS==level))
        differences <- diff(nonzeroIDX)
        kk <- which(differences>1 & differences<distanceThreshold)
        if(length(kk)){
          for (i in seq(1,length(kk))){
            DELETIONS[(nonzeroIDX[kk[i]]+1) : (nonzeroIDX[kk[i]+1]-1) ] <- level
          }
        }
      }}
    ## CLOSE GAPS IN BETWEEN TIGHT SVs DUPLICATIONS
    if(as.integer(max(DUPLICATIONS))>0){
      for(level in as.integer(max(DUPLICATIONS)):1){
        nonzeroIDX <- sort(which(DUPLICATIONS==level))
        differences <- diff(nonzeroIDX)
        kk <- which(differences>1 & differences<distanceThreshold)
        if(length(kk)){
          for (i in seq(1,length(kk))){
            DUPLICATIONS[(nonzeroIDX[kk[i]]+1) : (nonzeroIDX[kk[i]+1]-1) ] <- level
          }
        }
      }}
    ## CLOSE GAPS IN BETWEEN TIGHT SVs INSERTIONS
    if(as.integer(max(INSERTIONS))>0){
      for(level in as.integer(max(INSERTIONS)):1){
        nonzeroIDX <- sort(which(INSERTIONS==level))
        differences <- diff(nonzeroIDX)
        kk <- which(differences>1 & differences<distanceThreshold)
        if(length(kk)){
          for (i in seq(1,length(kk))){
            INSERTIONS[(nonzeroIDX[kk[i]]+1) : (nonzeroIDX[kk[i]+1]-1) ] <- level
          }
        }
      }}
    ## CLOSE GAPS IN BETWEEN TIGHT SVs INVERSIONS
    if(as.integer(max(INVERSIONS))>0){
      for(level in as.integer(max(INVERSIONS)):1){
        nonzeroIDX <- sort(which(INVERSIONS==level))
        differences <- diff(nonzeroIDX)
        kk <- which(differences>1 & differences<distanceThreshold)
        if(length(kk)){
          for (i in seq(1,length(kk))){
            INVERSIONS[(nonzeroIDX[kk[i]]+1) : (nonzeroIDX[kk[i]+1]-1) ] <- level
          }
        }
      }}
    
    #######################################################################################################################
    ## DELETIONS 
    #GET START and STOP of EVENTS - https://stackoverflow.com/questions/29184297/finding-the-start-and-stop-indices-in-sequence-in-r
    DELETIONS_DF <- data.frame(CHROM=character(),SVTYPE=character(),START=integer(), STOP=integer(),
                               LEN=integer(),NumOfCallers=integer(),Callers=character(),
                               cnproscanPortion=numeric(),delly2Portion=numeric(),lumpyPortion=numeric(),
                               cnvnatorPortion=numeric(),pindelPortion=numeric(),
                               cnproscanSubevents=integer(),delly2Subevents=integer(),lumpySubevents=integer(),
                               cnvnatorSubevents=integer(),pindelSubevents=integer()
    )
    if(as.integer(max(DELETIONS))>0){
      for(i in as.integer(max(DELETIONS)):1){
        # boundaries
        y <- which(DELETIONS>=i)
        START <- y[!(y-1) %in% y]
        STOP <- y[!(y+1) %in% y]
        #dataframe
        DELETIONS_DF_temp <- as.data.frame(cbind(START, STOP))
        DELETIONS_DF_temp$LEN <- DELETIONS_DF_temp$STOP - DELETIONS_DF_temp$START+1
        DELETIONS_DF_temp$NumOfCallers <- i
        DELETIONS_DF_temp$Callers <- ""
        DELETIONS_DF_temp$cnproscanPortion <- 0
        DELETIONS_DF_temp$delly2Portion <- 0
        DELETIONS_DF_temp$lumpyPortion <- 0
        DELETIONS_DF_temp$cnvnatorPortion <- 0
        DELETIONS_DF_temp$pindelPortion <- 0
        DELETIONS_DF_temp$cnproscanSubevents <- 0
        DELETIONS_DF_temp$delly2Subevents <- 0
        DELETIONS_DF_temp$lumpySubevents <- 0
        DELETIONS_DF_temp$cnvnatorSubevents <- 0
        DELETIONS_DF_temp$pindelSubevents <- 0
        # append 
        DELETIONS_DF <- rbind(DELETIONS_DF,DELETIONS_DF_temp)
      }
      
      # GET INFO
      for (j in 1:nrow(DELETIONS_DF)){
        # which callers called the region? and by how much
        callers=""
        if(any(which(cnproscanVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"cnproscan",sep = ",") # append to callers string
          DELETIONS_DF$cnproscanPortion[j] <- as.integer(length(which(cnproscanVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0))/(DELETIONS_DF$LEN[j]/100))
          tmp <- which(cnproscanVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]>0)
          DELETIONS_DF$cnproscanSubevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        if(any(which(dellyVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"delly2",sep = ",") # append to callers string
          DELETIONS_DF$delly2Portion[j] <- as.integer(length(which(dellyVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0))/(DELETIONS_DF$LEN[j]/100))
          tmp <- which(dellyVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0)
          DELETIONS_DF$delly2Subevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
          
        }
        if(any(which(lumpyVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"lumpy",sep = ",") # append to callers string
          DELETIONS_DF$lumpyPortion[j] <- as.integer(length(which(lumpyVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0))/(DELETIONS_DF$LEN[j]/100))
          tmp <- which(lumpyVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0)
          DELETIONS_DF$lumpySubevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        if(any(which(cnvnatorVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"cnvnator",sep = ",") # append to callers string
          DELETIONS_DF$cnvnatorPortion[j] <- as.integer(length(which(cnvnatorVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0))/(DELETIONS_DF$LEN[j]/100))
          tmp <- which(cnvnatorVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0)
          DELETIONS_DF$cnvnatorSubevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        if(any(which(pindelVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"pindel",sep = ",") # append to callers string
          DELETIONS_DF$pindelPortion[j] <- as.integer(length(which(pindelVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0))/(DELETIONS_DF$LEN[j]/100))
          tmp <- which(pindelVEC_DEL[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])]!=0)
          DELETIONS_DF$pindelSubevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        callers <- gsub("^,","",callers)
        DELETIONS_DF$Callers[j] <- callers
      }
      
      DELETIONS_DF$SVTYPE <- "DEL"
      DELETIONS_DF$CHROM <- refHeader
    }
    
    DELETIONS_DF <- DELETIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                    "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                    "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )]
    
    colnames(DELETIONS_DF) <- c("CHROMOSOME","START","STOP","LENGTH_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )
    
    #######################################################################################################################
    ## DUPLICATIONS
    DUPLICATIONS_DF <- data.frame(CHROM=character(),SVTYPE=character(),START=integer(), STOP=integer(),
                                  LEN=integer(),NumOfCallers=integer(),Callers=character(),
                                  cnproscanPortion=numeric(),delly2Portion=numeric(),lumpyPortion=numeric(),
                                  cnvnatorPortion=numeric(),pindelPortion=numeric(),
                                  cnproscanSubevents=integer(),delly2Subevents=integer(),lumpySubevents=integer(),
                                  cnvnatorSubevents=integer(),pindelSubevents=integer()
    )
    if(as.integer(max(DUPLICATIONS))>0){
      for(i in as.integer(max(DUPLICATIONS)):1){
        # boundaries
        y <- which(DUPLICATIONS>=i)
        START <- y[!(y-1) %in% y]
        STOP <- y[!(y+1) %in% y]
        #dataframe
        DUPLICATIONS_DF_temp <- as.data.frame(cbind(START, STOP))
        DUPLICATIONS_DF_temp$LEN <- DUPLICATIONS_DF_temp$STOP - DUPLICATIONS_DF_temp$START+1
        DUPLICATIONS_DF_temp$NumOfCallers <- i
        DUPLICATIONS_DF_temp$Callers <- ""
        DUPLICATIONS_DF_temp$cnproscanPortion <- 0
        DUPLICATIONS_DF_temp$delly2Portion <- 0
        DUPLICATIONS_DF_temp$lumpyPortion <- 0
        DUPLICATIONS_DF_temp$cnvnatorPortion <- 0
        DUPLICATIONS_DF_temp$pindelPortion <- 0
        DUPLICATIONS_DF_temp$cnproscanSubevents <- 0
        DUPLICATIONS_DF_temp$delly2Subevents <- 0
        DUPLICATIONS_DF_temp$lumpySubevents <- 0
        DUPLICATIONS_DF_temp$cnvnatorSubevents <- 0
        DUPLICATIONS_DF_temp$pindelSubevents <- 0
        # append 
        DUPLICATIONS_DF <- rbind(DUPLICATIONS_DF,DUPLICATIONS_DF_temp)
      }
      
      # GET INFO
      for (j in 1:nrow(DUPLICATIONS_DF)){
        # which callers called the region? and by how much
        callers=""
        if(any(which(cnproscanVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"cnproscan",sep = ",") # append to callers string
          DUPLICATIONS_DF$cnproscanPortion[j] <- as.integer(length(which(cnproscanVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0))/(DUPLICATIONS_DF$LEN[j]/100))
          tmp <- which(cnproscanVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]>0)
          DUPLICATIONS_DF$cnproscanSubevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        if(any(which(dellyVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"delly2",sep = ",") # append to callers string
          DUPLICATIONS_DF$delly2Portion[j] <- as.integer(length(which(dellyVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0))/(DUPLICATIONS_DF$LEN[j]/100))
          tmp <- which(dellyVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0)
          DUPLICATIONS_DF$delly2Subevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
          
        }
        if(any(which(lumpyVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"lumpy",sep = ",") # append to callers string
          DUPLICATIONS_DF$lumpyPortion[j] <- as.integer(length(which(lumpyVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0))/(DUPLICATIONS_DF$LEN[j]/100))
          tmp <- which(lumpyVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0)
          DUPLICATIONS_DF$lumpySubevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        if(any(which(cnvnatorVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"cnvnator",sep = ",") # append to callers string
          DUPLICATIONS_DF$cnvnatorPortion[j] <- as.integer(length(which(cnvnatorVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0))/(DUPLICATIONS_DF$LEN[j]/100))
          tmp <- which(cnvnatorVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0)
          DUPLICATIONS_DF$cnvnatorSubevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        if(any(which(pindelVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"pindel",sep = ",") # append to callers string
          DUPLICATIONS_DF$pindelPortion[j] <- as.integer(length(which(pindelVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0))/(DUPLICATIONS_DF$LEN[j]/100))
          tmp <- which(pindelVEC_DUP[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])]!=0)
          DUPLICATIONS_DF$pindelSubevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        callers <- gsub("^,","",callers)
        DUPLICATIONS_DF$Callers[j] <- callers
      }
      DUPLICATIONS_DF$SVTYPE <- "DUP"
      DUPLICATIONS_DF$CHROM <- refHeader
    }
    DUPLICATIONS_DF <- DUPLICATIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                          "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                          "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )]
    colnames(DUPLICATIONS_DF) <- c("CHROMOSOME","START","STOP","LENGTH_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                   "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                   "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents"  )
    
    #######################################################################################################################
    ## INSERTIONS
    INSERTIONS_DF <- data.frame(CHROM=character(),SVTYPE=character(),START=integer(), STOP=integer(),
                                LEN=integer(),NumOfCallers=integer(),Callers=character(),
                                cnproscanPortion=numeric(),delly2Portion=numeric(),lumpyPortion=numeric(),
                                cnvnatorPortion=numeric(),pindelPortion=numeric(),
                                cnproscanSubevents=integer(),delly2Subevents=integer(),lumpySubevents=integer(),
                                cnvnatorSubevents=integer(),pindelSubevents=integer()
    )
    if(as.integer(max(INSERTIONS))>0){
      for(i in as.integer(max(INSERTIONS)):1){
        # boundaries
        y <- which(INSERTIONS>=i)
        START <- y[!(y-1) %in% y]
        STOP <- y[!(y+1) %in% y]
        #dataframe
        INSERTIONS_DF_temp <- as.data.frame(cbind(START, STOP))
        INSERTIONS_DF_temp$LEN <- INSERTIONS_DF_temp$STOP - INSERTIONS_DF_temp$START+1
        INSERTIONS_DF_temp$NumOfCallers <- i
        INSERTIONS_DF_temp$Callers <- ""
        INSERTIONS_DF_temp$cnproscanPortion <- 0
        INSERTIONS_DF_temp$delly2Portion <- 0
        INSERTIONS_DF_temp$lumpyPortion <- 0
        INSERTIONS_DF_temp$cnvnatorPortion <- 0
        INSERTIONS_DF_temp$pindelPortion <- 0
        INSERTIONS_DF_temp$cnproscanSubevents <- 0
        INSERTIONS_DF_temp$delly2Subevents <- 0
        INSERTIONS_DF_temp$lumpySubevents <- 0
        INSERTIONS_DF_temp$cnvnatorSubevents <- 0
        INSERTIONS_DF_temp$pindelSubevents <- 0
        # append 
        INSERTIONS_DF <- rbind(INSERTIONS_DF,INSERTIONS_DF_temp)
      }
      
      # GET INFO
      for (j in 1:nrow(INSERTIONS_DF)){
        # which callers called the region? and by how much
        callers=""
        
        if(any(which(dellyVEC_INS[(INSERTIONS_DF$START[j]):(INSERTIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"delly2",sep = ",") # append to callers string
          INSERTIONS_DF$delly2Portion[j] <- as.integer(length(which(dellyVEC_INS[(INSERTIONS_DF$START[j]):(INSERTIONS_DF$STOP[j])]!=0))/(INSERTIONS_DF$LEN[j]/100))
          tmp <- which(dellyVEC_INS[(INSERTIONS_DF$START[j]):(INSERTIONS_DF$STOP[j])]!=0)
          INSERTIONS_DF$delly2Subevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        if(any(which(pindelVEC_INS[(INSERTIONS_DF$START[j]):(INSERTIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"pindel",sep = ",") # append to callers string
          INSERTIONS_DF$pindelPortion[j] <- as.integer(length(which(pindelVEC_INS[(INSERTIONS_DF$START[j]):(INSERTIONS_DF$STOP[j])]!=0))/(INSERTIONS_DF$LEN[j]/100))
          tmp <- which(pindelVEC_INS[(INSERTIONS_DF$START[j]):(INSERTIONS_DF$STOP[j])]!=0)
          INSERTIONS_DF$pindelSubevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        callers <- gsub("^,","",callers)
        INSERTIONS_DF$Callers[j] <- callers
      }
      INSERTIONS_DF$SVTYPE <- "INS"
      INSERTIONS_DF$CHROM <- refHeader
    }
    INSERTIONS_DF <- INSERTIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                      "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                      "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )]
    colnames(INSERTIONS_DF) <- c("CHROMOSOME","START","STOP","LENGTH_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                 "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                 "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )
    
    #######################################################################################################################
    ##INVERSIONS
    INVERSIONS_DF <- data.frame(CHROM=character(),SVTYPE=character(),START=integer(), STOP=integer(),
                                LEN=integer(),NumOfCallers=integer(),Callers=character(),
                                cnproscanPortion=numeric(),delly2Portion=numeric(),lumpyPortion=numeric(),
                                cnvnatorPortion=numeric(),pindelPortion=numeric(),
                                cnproscanSubevents=integer(),delly2Subevents=integer(),lumpySubevents=integer(),
                                cnvnatorSubevents=integer(),pindelSubevents=integer()
    )
    if(as.integer(max(INVERSIONS))>0){
      for(i in as.integer(max(INVERSIONS)):1){
        # boundaries
        y <- which(INVERSIONS>=i)
        START <- y[!(y-1) %in% y]
        STOP <- y[!(y+1) %in% y]
        #dataframe
        INVERSIONS_DF_temp <- as.data.frame(cbind(START, STOP))
        INVERSIONS_DF_temp$LEN <- INVERSIONS_DF_temp$STOP - INVERSIONS_DF_temp$START+1
        INVERSIONS_DF_temp$NumOfCallers <- i
        INVERSIONS_DF_temp$Callers <- ""
        INVERSIONS_DF_temp$cnproscanPortion <- 0
        INVERSIONS_DF_temp$delly2Portion <- 0
        INVERSIONS_DF_temp$lumpyPortion <- 0
        INVERSIONS_DF_temp$cnvnatorPortion <- 0
        INVERSIONS_DF_temp$pindelPortion <- 0
        INVERSIONS_DF_temp$cnproscanSubevents <- 0
        INVERSIONS_DF_temp$delly2Subevents <- 0
        INVERSIONS_DF_temp$lumpySubevents <- 0
        INVERSIONS_DF_temp$cnvnatorSubevents <- 0
        INVERSIONS_DF_temp$pindelSubevents <- 0
        # append 
        INVERSIONS_DF <- rbind(INVERSIONS_DF,INVERSIONS_DF_temp)
      }
      
      # GET INFO
      for (j in 1:nrow(INVERSIONS_DF)){
        # which callers called the region? and by how much
        callers=""
        
        if(any(which(dellyVEC_INV[(INVERSIONS_DF$START[j]):(INVERSIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"delly2",sep = ",") # append to callers string
          INVERSIONS_DF$delly2Portion[j] <- as.integer(length(which(dellyVEC_INV[(INVERSIONS_DF$START[j]):(INVERSIONS_DF$STOP[j])]!=0))/(INVERSIONS_DF$LEN[j]/100))
          tmp <- which(dellyVEC_INV[(INVERSIONS_DF$START[j]):(INVERSIONS_DF$STOP[j])]!=0)
          INVERSIONS_DF$delly2Subevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        if(any(which(pindelVEC_INV[(INVERSIONS_DF$START[j]):(INVERSIONS_DF$STOP[j])]!=0))){
          callers <- paste(callers,"pindel",sep = ",") # append to callers string
          INVERSIONS_DF$pindelPortion[j] <- as.integer(length(which(pindelVEC_INV[(INVERSIONS_DF$START[j]):(INVERSIONS_DF$STOP[j])]!=0))/(INVERSIONS_DF$LEN[j]/100))
          tmp <- which(pindelVEC_INV[(INVERSIONS_DF$START[j]):(INVERSIONS_DF$STOP[j])]!=0)
          INVERSIONS_DF$pindelSubevents[j] <- as.integer(length(tmp[!(tmp-1) %in% tmp]))
        }
        callers <- gsub("^,","",callers)
        INVERSIONS_DF$Callers[j] <- callers
      }
      INVERSIONS_DF$SVTYPE <- "INV"
      INVERSIONS_DF$CHROM <- refHeader
    }
    
    INVERSIONS_DF <- INVERSIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                      "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                      "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents")]
    colnames(INVERSIONS_DF) <- c("CHROMOSOME","START","STOP","LENGTH_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                 "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                 "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents")
    
    #######################################################################################################################
    # GET TOGETHER # WHAT ABOUT LARGE FILES?
    SHORT_RESULTS <- rbindlist(list(SHORT_RESULTS,DELETIONS_DF, DUPLICATIONS_DF, INSERTIONS_DF, INVERSIONS_DF),fill=TRUE)
    
  }
  
  #######################################################################################################################
  #sort
  SHORT_RESULTS <- SHORT_RESULTS[order(SHORT_RESULTS$CHROMOSOME, SHORT_RESULTS$START),]
  
  #write TABLE
  dir.create(file.path(getwd(),dirname(output)), recursive = TRUE)
  write.table(SHORT_RESULTS, file=output, sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE )
  
  
  #######################################################################################################################
  # VENN DIAGRAM of CALLERS
  SHORT_RESULTS$ID <- paste(SHORT_RESULTS$CHROMOSOME,SHORT_RESULTS$START,SHORT_RESULTS$STOP,SHORT_RESULTS$LENGTH_MERGED,SHORT_RESULTS$SVTYPE_MERGED,sep="_")
  
  
  ## detect files by caller
  delly_SVs <- SHORT_RESULTS[which(grepl( "delly2", SHORT_RESULTS$Callers, fixed = TRUE)),"ID"]
  lumpy_SVs <- SHORT_RESULTS[which(grepl( "lumpy",  SHORT_RESULTS$Callers, fixed = TRUE)),"ID"]
  cnvnator_SVs <- SHORT_RESULTS[which(grepl( "cnvnator",  SHORT_RESULTS$Callers, fixed = TRUE)),"ID"]
  pindel_SVs <- SHORT_RESULTS[which(grepl( "pindel",  SHORT_RESULTS$Callers, fixed = TRUE)),"ID"]
  cnproscan_SVs <- SHORT_RESULTS[which(grepl( "cnproscan",  SHORT_RESULTS$Callers, fixed = TRUE)),"ID"]
  
  x <- list(
    DELLY = unique(delly_SVs$ID),
    LUMPY = unique(lumpy_SVs$ID),
    CNVNATOR = unique(cnvnator_SVs$ID),
    PINDEL = unique(pindel_SVs$ID),
    CNPROSCAN = unique(cnproscan_SVs$ID)
  )
  venn <- Venn(x)
  data <- process_data(venn)
  
  png(file=venn_png,width=30,height=20,units="cm",res=300)
  p <- ggplot() +
    # 1. region count layer
    geom_sf(aes(fill = count), data = venn_region(data)) +
    # 2. set edge layer
    geom_sf(color=c("blue","red","yellow","green","purple"), alpha = 0.7, data = venn_setedge(data), show.legend = TRUE) +
    # 3. set label layer
    geom_sf_text(aes(label = c("DELLY","LUMPY","CNVNATOR","PINDEL","CNPROSCAN")), data = venn_setlabel(data)) +
    # 4. region label layer
    geom_sf_label(aes(label = count), color="white", label.size  = NA, size=5, alpha = 0.0, data = venn_region(data)) +
    # 5. color theme
    scale_color_brewer(palette = "Set2") +
    labs(title = "Venn diagram of detected SVs by callers") + theme_void()
  print(p)
  dev.off()
  Sys.sleep(5)
  
  #######################################################################################################################
  # PIE PLOT of SV TYPES
  SHORT_RESULTS$DELLY <- 0
  SHORT_RESULTS$LUMPY <- 0
  SHORT_RESULTS$PINDEL <- 0
  SHORT_RESULTS$CNVNATOR <- 0 
  SHORT_RESULTS$CNPROSCAN <- 0 
  
  SHORT_RESULTS[which(grepl( "delly2", SHORT_RESULTS$Callers, fixed = TRUE)),"DELLY"] <- 1
  SHORT_RESULTS[which(grepl( "lumpy",  SHORT_RESULTS$Callers, fixed = TRUE)),"LUMPY"] <- 1
  SHORT_RESULTS[which(grepl( "pindel",  SHORT_RESULTS$Callers, fixed = TRUE)),"PINDEL"] <- 1
  SHORT_RESULTS[which(grepl( "cnvnator",  SHORT_RESULTS$Callers, fixed = TRUE)),"CNVNATOR"] <- 1
  SHORT_RESULTS[which(grepl( "cnproscan",  SHORT_RESULTS$Callers, fixed = TRUE)),"CNPROSCAN"]   <- 1
  
  data <- data.frame(
    type=c(rep("DUP",5), rep("DEL",5), rep("INV",2), rep("INS",2)      ),
    caller= c("DELLY","LUMPY","PINDEL","CNVNATOR","CNPROSCAN",
              "DELLY","LUMPY","PINDEL","CNVNATOR","CNPROSCAN",
              "DELLY","PINDEL",
              "DELLY","PINDEL"),
    value=rep(0,times=14)
  )
  
  #DUP
  data$value[1] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="DUP" & SHORT_RESULTS$DELLY==1 )
  data$value[2] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="DUP" & SHORT_RESULTS$LUMPY==1 )
  data$value[3] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="DUP" & SHORT_RESULTS$PINDEL==1 )
  data$value[4] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="DUP" & SHORT_RESULTS$CNVNATOR==1 )
  data$value[5] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="DUP" & SHORT_RESULTS$CNPROSCAN==1 )
  #DEL
  data$value[6] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="DEL" & SHORT_RESULTS$DELLY==1 )
  data$value[7] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="DEL" & SHORT_RESULTS$LUMPY==1 )
  data$value[8] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="DEL" & SHORT_RESULTS$PINDEL==1 )
  data$value[9] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="DEL" & SHORT_RESULTS$CNVNATOR==1 )
  data$value[10] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="DEL" & SHORT_RESULTS$CNPROSCAN==1 )
  #INV
  data$value[11] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="INV" & SHORT_RESULTS$DELLY==1 )
  data$value[12] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="INV" & SHORT_RESULTS$PINDEL==1 )
  #INS
  data$value[13] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="INS" & SHORT_RESULTS$DELLY==1 )
  data$value[14] <- sum(SHORT_RESULTS$SVTYPE_MERGED=="INS" & SHORT_RESULTS$PINDEL==1 )
  
  png(file=barchart_png,width=40,height=20,units="cm",res=300)
  p <- ggplot(data, aes(x = type, y = value, fill = caller)) +   geom_bar(stat = "identity") + 
    xlab("SV type") + ylab("Count") + labs(title = ("Barchart of represented SV types"))
  print(p)
  dev.off()
  Sys.sleep(5)
  
}





# develop and test
# args <- c("resultsTEST/vcf_merged/SRR20631331.procaryaSV_merge0.tsv",
#           "results/references/Vibrio_cholerae_strain_A19.fasta",
#           NA,NA,NA,
#           "results/lumpy/SRR20631331/SRR20631331.vcf","results/delly2/SRR20631331/SRR20631331.vcf",
#           "results/cnvnator/SRR20631331/SRR20631331.vcf","results/pindel/SRR20631331/SRR20631331.vcf",
#           "results/cnproscan/SRR20631331/SRR20631331.vcf")
# setwd("/home/rj/1TB/ProcaryaSV_test_multichrom/")

#run as Rscript

args <- commandArgs(trailingOnly = T)
run_all(args)