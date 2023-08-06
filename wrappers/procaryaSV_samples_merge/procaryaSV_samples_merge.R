
library(data.table)
library(IRanges)
library(stringr)

## test
setwd("/home/rj/4TB/PHD_TESTING/K.pneumonaie2/")
path <- "/home/rj/4TB/PHD_TESTING/K.pneumonaie2/results/merged_procaryaSV"
tsv.names <- dir(path, pattern =".tsv$",recursive = TRUE)

input.files <- paste(path,tsv.names,sep = "/")

allCNVs <- rbindlist(lapply(input.files,fread))

output <- "testOUTPUT.tsv"
heatmap_png <-"testheatmap.png"
fastaFile <- "results/references/KP_ref.fasta"
min_sv_length <- 100
max_sv_length <- NA
distanceThreshold <- 2000
minCallers <- 2
minSamples <- 2


run_all <- function(args){
  # arguments
  output <- args[1]
  heatmap_png <- args[2]
  fastaFile <- args[3] 
  min_sv_length <- as.integer(args[4]) #default 1 
  max_sv_length <- as.integer(args[5]) #default refLength/3
  distanceThreshold <- as.integer(args[6]) 
  minCallers <- as.integer(args[7]) 
  minSamples <- as.integer(args[9]) 
  input.files <- args[10:length(args)]
  
  ############################################################################################################################################################
  ## DEFAULT VALUES
  if(is.na(distanceThreshold)){distanceThreshold <- 2000} #default
  if(is.na(min_sv_length)){min_sv_length <- 50} #default
  # max_sv_length depends on chromosome length
  
  ## read FASTA
  referenceLengthList <- seqinr::getLength(seqinr::read.fasta(fastaFile, seqonly=TRUE))
  numberOfFastaContigs <- length(referenceLengthList)
  fastaContiglist <- seqinr::getName(seqinr::read.fasta(fastaFile, seqonly=FALSE))
  
  threshold <- round(length(input.files)/10)#10%
  if(threshold==0){threshold <- 1}
  ##############################################################################
  ## LOAD TSV FILES
  allCNVs <- rbindlist(lapply(input.files,fread))
  allCNVs <- allCNVs[allCNVs$MinSupport>=minCallers,]
  
  ############################################################################################################################################################
  ############################################################################################################################################################
  ## GO OVER CHROMOSOMES
  
  SHORT_RESULTS <- data.table()
  
  contigNum <- 1
  
  for (contigNum in 1:numberOfFastaContigs){
    
    refHeader <- fastaContiglist[contigNum]
    refLength <- referenceLengthList[contigNum]
    if(is.na(max_sv_length)){max_sv_length <- refLength/3} #default
    
    # GET CURRENT CHROMOSOME and SVs that fit
    TMPallCNVs <- allCNVs[allCNVs$CHROMOSOME==refHeader,]
    TMPallCNVs <- allCNVs[allCNVs$LENGTH_MERGED<=max_sv_length,]  #
    TMPallCNVs <- allCNVs[allCNVs$LENGTH_MERGED>=min_sv_length,]
    
    DELETIONS <- rep(0,refLength)
    DUPLICATIONS <- rep(0,refLength)
    INSERTIONS <- rep(0,refLength)
    INVERSIONS <- rep(0,refLength)
    
    # SVs SIGNAL
    if(nrow(TMPallCNVs)){
      for (j in 1:nrow(TMPallCNVs)){
        if(TMPallCNVs$SVTYPE[j]=="DUP"){DUPLICATIONS[TMPallCNVs$START[j]:TMPallCNVs$STOP[j]] <- DUPLICATIONS[TMPallCNVs$START[j]:TMPallCNVs$STOP[j]]+1}
        if(TMPallCNVs$SVTYPE[j]=="DEL"){DELETIONS[TMPallCNVs$START[j]:TMPallCNVs$STOP[j]] <- DELETIONS[TMPallCNVs$START[j]:TMPallCNVs$STOP[j]]+1}
        if(TMPallCNVs$SVTYPE[j]=="INS"){INSERTIONS[TMPallCNVs$START[j]:TMPallCNVs$STOP[j]] <- INSERTIONS[TMPallCNVs$START[j]:TMPallCNVs$STOP[j]]+1}
        if(TMPallCNVs$SVTYPE[j]=="INV"){INVERSIONS[TMPallCNVs$START[j]:TMPallCNVs$STOP[j]] <- INVERSIONS[TMPallCNVs$START[j]:TMPallCNVs$STOP[j]]+1}
      }}
    
    
    DUPLICATIONS[is.na(DUPLICATIONS)] <- 0
    DELETIONS[is.na(DELETIONS)] <- 0
    INSERTIONS[is.na(INSERTIONS)] <- 0
    INVERSIONS[is.na(INVERSIONS)] <- 0
    
    #######################################################################################################################
    ## DELETIONS  SAMPLES MERGE
    DELETIONS_DF <- data.frame(ID=character(),CHROMOSOME=character(),START=integer(), STOP=integer(),
                               LENGTH_MERGED=integer(),SVTYPE_MERGED=character(), MinSupport=integer(),MaxSupport=integer(),
                               maxSupSTART=integer(),maxSupSTOP=integer(),
                               Samples=character()
    )
    DELETIONS_DF_tmp <- DELETIONS_DF
    
    
    if(as.integer(max(DELETIONS))>0){
      
      #A get regions boundaries
      for(level in seq(as.integer(max(DELETIONS)),threshold,by=-1) ){
        # boundaries https://stackoverflow.com/questions/29184297/finding-the-start-and-stop-indices-in-sequence-in-r
        y <- which(DELETIONS>=level)
        START <- y[!(y-1) %in% y]
        STOP <- y[!(y+1) %in% y]
        # dataframe filling
        ID <- c(1:length(STOP)) 
        DF_temp <- as.data.frame(cbind(ID,START, STOP))
        DF_temp$LENGTH_MERGED <- DF_temp$STOP - DF_temp$START+1
        DF_temp$MinSupport <- level
        DF_temp$MaxSupport <- 0
        DF_temp$maxSupSTART <- 0
        DF_temp$maxSupSTOP <- 0
        DF_temp$Samples <- ""
        # append 
        DELETIONS_DF <- rbind(DELETIONS_DF,DF_temp)
      }
      
      # modify DF columns
      DELETIONS_DF$CHROMOSOME <- refHeader
      DELETIONS_DF$ID <- c(1:nrow(DELETIONS_DF))
      DELETIONS_DF$ID <- paste(contigNum,"DEL",DELETIONS_DF$ID,sep="_")
      DELETIONS_DF$SVTYPE_MERGED <- "DEL"
      
    }
    
    #B find overlapping events
    if(nrow(DELETIONS_DF)>0){
      
      for (j in 1:nrow(DELETIONS_DF)){
        # https://stackoverflow.com/questions/37754509/finding-overlapping-intervals
        DUPLICATES <- DELETIONS_DF[subjectHits(IRanges::findOverlaps(IRanges(DELETIONS_DF$START[j], DELETIONS_DF$STOP[j]), IRanges(DELETIONS_DF$START, DELETIONS_DF$STOP), type="equal", maxgap = distanceThreshold) ), "ID"]
        DELETIONS_DF$Samples[j] <- as.character(paste(DUPLICATES, collapse = ";"))
      }
      rm(DUPLICATES)
      
      #B collapse overlapping events
      used <- character()
      for (j in 1:nrow(DELETIONS_DF)){
        if(!DELETIONS_DF$ID[j] %in% used){
          # get matches
          matches <- unlist(str_split(DELETIONS_DF$Samples[j],";"))
          
          if(length(matches)!=0){
            # df
            temp <- DELETIONS_DF[which(DELETIONS_DF$ID %in% matches),]
            # find max support CNV and use it
            idx <- which.max(DELETIONS_DF$MinSupport)
            temp$maxSupSTART[idx] <- temp$START[idx]
            temp$maxSupSTOP[idx] <- temp$STOP[idx]
            temp$MaxSupport[idx] <- max(temp$MinSupport)
            temp$MinSupport[idx] <- min(temp$MinSupport)
            temp$START[idx] <- min(temp$START)
            temp$STOP[idx] <- max(temp$STOP)
            #delete from pool
            used <- append(used, temp$ID)
            temp <- temp[idx,]
            
            #append new
            DELETIONS_DF_tmp <- rbind(DELETIONS_DF_tmp,temp)
            
          }else(
            #append current
            DELETIONS_DF_tmp <- rbind(DELETIONS_DF_tmp, DELETIONS_DF[j,])
          )
          
        }
      }
      
      #replace and get new IDs
      DELETIONS_DF <- DELETIONS_DF_tmp
      DELETIONS_DF$ID <- c(1:nrow(DELETIONS_DF))
      DELETIONS_DF$ID <- paste(contigNum,"DEL",DELETIONS_DF$ID,sep="_")
      DELETIONS_DF$Samples <- ""
      
      
      #C BACKTRACKING - FIND IN WHICH SAMPLES THEY ARE
      TMPx <- TMPallCNVs[TMPallCNVs$SVTYPE_MERGED=="DEL",]
      DELETIONS_DF$NumSamples <- 0
      # DELETIONS_DF$IDs <- ""
      for (j in 1:nrow(DELETIONS_DF)){
        # which samples overlap with the region?
        MATCHES <- TMPx[subjectHits(IRanges::findOverlaps(IRanges(DELETIONS_DF$START[j], DELETIONS_DF$STOP[j]), IRanges(TMPx$START, TMPx$STOP),type="equal",maxgap = distanceThreshold) ),]
        DELETIONS_DF$Samples[j] <- paste(as.character(unique(MATCHES$SAMPLE)), collapse = ";")
        DELETIONS_DF$NumSamples[j] <- length(unique(MATCHES$SAMPLE))
        # DELETIONS_DF$IDs[j] <- paste(as.character(unique(MATCHES$ID)), collapse = ";")
      }
      
      DELETIONS_DF <- DELETIONS_DF[,c("ID","CHROMOSOME","START","STOP","LENGTH_MERGED","SVTYPE_MERGED", "NumSamples", "Samples")]
      DELETIONS_DF <- DELETIONS_DF[DELETIONS_DF$NumSamples>=minSamples,]
      
    }
    
    
    
    
    
    
    #######################################################################################################################
    ## DUPLICATIONS
    DUPLICATIONS_DF <- data.frame(ID=character(),CHROMOSOME=character(),START=integer(), STOP=integer(),
                                  LENGTH_MERGED=integer(),SVTYPE_MERGED=character(), MinSupport=integer(),MaxSupport=integer(),
                                  maxSupSTART=integer(),maxSupSTOP=integer(),
                                  Samples=character()
    )
    DUPLICATIONS_DF_tmp <- DUPLICATIONS_DF
    
    
    if(as.integer(max(DUPLICATIONS))>0){
      
      #A get regions boundaries
      for(level in seq(as.integer(max(DUPLICATIONS)),threshold,by=-1) ){
        # boundaries https://stackoverflow.com/questions/29184297/finding-the-start-and-stop-indices-in-sequence-in-r
        y <- which(DUPLICATIONS>=level)
        START <- y[!(y-1) %in% y]
        STOP <- y[!(y+1) %in% y]
        # dataframe filling
        ID <- c(1:length(STOP)) 
        DF_temp <- as.data.frame(cbind(ID,START, STOP))
        DF_temp$LENGTH_MERGED <- DF_temp$STOP - DF_temp$START+1
        DF_temp$MinSupport <- level
        DF_temp$MaxSupport <- 0
        DF_temp$maxSupSTART <- 0
        DF_temp$maxSupSTOP <- 0
        DF_temp$Samples <- ""
        # append 
        DUPLICATIONS_DF <- rbind(DUPLICATIONS_DF,DF_temp)
      }
      
      # modify DF columns
      DUPLICATIONS_DF$CHROMOSOME <- refHeader
      DUPLICATIONS_DF$ID <- c(1:nrow(DUPLICATIONS_DF))
      DUPLICATIONS_DF$ID <- paste(contigNum,"DUP",DUPLICATIONS_DF$ID,sep="_")
      DUPLICATIONS_DF$SVTYPE_MERGED <- "DUP"
    }
    
    #B find overlapping events
    if(nrow(DUPLICATIONS_DF)>0){
      
      for (j in 1:nrow(DUPLICATIONS_DF)){
        # https://stackoverflow.com/questions/37754509/finding-overlapping-intervals
        DUPLICATES <- DUPLICATIONS_DF[subjectHits(IRanges::findOverlaps(IRanges(DUPLICATIONS_DF$START[j], DUPLICATIONS_DF$STOP[j]), IRanges(DUPLICATIONS_DF$START, DUPLICATIONS_DF$STOP), type="equal", maxgap = distanceThreshold) ), "ID"]
        DUPLICATIONS_DF$Samples[j] <- as.character(paste(DUPLICATES, collapse = ";"))
      }
      rm(DUPLICATES)
      
      #B collapse overlapping events
      used <- character()
      for (j in 1:nrow(DUPLICATIONS_DF)){
        if(!DUPLICATIONS_DF$ID[j] %in% used){
          # get matches
          matches <- unlist(str_split(DUPLICATIONS_DF$Samples[j],";"))
          
          if(length(matches)!=0){
            # df
            temp <- DUPLICATIONS_DF[which(DUPLICATIONS_DF$ID %in% matches),]
            # find max support CNV and use it
            idx <- which.max(DUPLICATIONS_DF$MinSupport)
            temp$maxSupSTART[idx] <- temp$START[idx]
            temp$maxSupSTOP[idx] <- temp$STOP[idx]
            temp$MaxSupport[idx] <- max(temp$MinSupport)
            temp$MinSupport[idx] <- min(temp$MinSupport)
            temp$START[idx] <- min(temp$START)
            temp$STOP[idx] <- max(temp$STOP)
            #delete from pool
            used <- append(used, temp$ID)
            temp <- temp[idx,]
            
            #append new
            DUPLICATIONS_DF_tmp <- rbind(DUPLICATIONS_DF_tmp,temp)
            
          }else(
            #append current
            DUPLICATIONS_DF_tmp <- rbind(DUPLICATIONS_DF_tmp, DUPLICATIONS_DF[j,])
          )
          
        }
      }
      
      #replace and get new IDs
      DUPLICATIONS_DF <- DUPLICATIONS_DF_tmp
      DUPLICATIONS_DF$ID <- c(1:nrow(DUPLICATIONS_DF))
      DUPLICATIONS_DF$ID <- paste(contigNum,"DUP",DUPLICATIONS_DF$ID,sep="_")
      DUPLICATIONS_DF$Samples <- ""
      
      
      #C BACKTRACKING - FIND IN WHICH SAMPLES THEY ARE
      TMPx <- TMPallCNVs[TMPallCNVs$SVTYPE_MERGED=="DUP",]
      DUPLICATIONS_DF$NumSamples <- 0
      # DUPLICATIONS_DF$IDs <- ""
      for (j in 1:nrow(DUPLICATIONS_DF)){
        # which samples overlap with the region?
        MATCHES <- TMPx[subjectHits(IRanges::findOverlaps(IRanges(DUPLICATIONS_DF$START[j], DUPLICATIONS_DF$STOP[j]), IRanges(TMPx$START, TMPx$STOP),type="equal",maxgap = distanceThreshold) ),]
        DUPLICATIONS_DF$Samples[j] <- paste(as.character(unique(MATCHES$SAMPLE)), collapse = ";")
        DUPLICATIONS_DF$NumSamples[j] <- length(unique(MATCHES$SAMPLE))
        # DUPLICATIONS_DF$IDs[j] <- paste(as.character(unique(MATCHES$ID)), collapse = ";")
      }
      DUPLICATIONS_DF <- DUPLICATIONS_DF[,c("ID","CHROMOSOME","START","STOP","LENGTH_MERGED","SVTYPE_MERGED", "NumSamples", "Samples")]
      DUPLICATIONS_DF <- DUPLICATIONS_DF[DUPLICATIONS_DF$NumSamples>=minSamples,]
    }
    
    
    
    
    
    
    
    #######################################################################################################################
    ## INSERTIONS
    INSERTIONS_DF <- data.frame(ID=character(),CHROMOSOME=character(),START=integer(), STOP=integer(),
                                LENGTH_MERGED=integer(),SVTYPE_MERGED=character(), MinSupport=integer(),MaxSupport=integer(),
                                maxSupSTART=integer(),maxSupSTOP=integer(),
                                Samples=character()
    )
    INSERTIONS_DF_tmp <- INSERTIONS_DF
    
    
    if(as.integer(max(INSERTIONS))>0){
      
      #A get regions boundaries
      for(level in seq(as.integer(max(INSERTIONS)),threshold,by=-1) ){
        # boundaries https://stackoverflow.com/questions/29184297/finding-the-start-and-stop-indices-in-sequence-in-r
        y <- which(INSERTIONS>=level)
        START <- y[!(y-1) %in% y]
        STOP <- y[!(y+1) %in% y]
        # dataframe filling
        ID <- c(1:length(STOP)) 
        DF_temp <- as.data.frame(cbind(ID,START, STOP))
        DF_temp$LENGTH_MERGED <- DF_temp$STOP - DF_temp$START+1
        DF_temp$MinSupport <- level
        DF_temp$MaxSupport <- 0
        DF_temp$maxSupSTART <- 0
        DF_temp$maxSupSTOP <- 0
        DF_temp$Samples <- ""
        # append 
        INSERTIONS_DF <- rbind(INSERTIONS_DF,DF_temp)
      }
      
      # modify DF columns
      INSERTIONS_DF$CHROMOSOME <- refHeader
      INSERTIONS_DF$ID <- c(1:nrow(INSERTIONS_DF))
      INSERTIONS_DF$ID <- paste(contigNum,"INS",INSERTIONS_DF$ID,sep="_")
      INSERTIONS_DF$SVTYPE_MERGED <- "INS"
      
    }
    
    #B find overlapping events
    if(nrow(INSERTIONS_DF)>0){
      
      for (j in 1:nrow(INSERTIONS_DF)){
        # https://stackoverflow.com/questions/37754509/finding-overlapping-intervals
        DUPLICATES <- INSERTIONS_DF[subjectHits(IRanges::findOverlaps(IRanges(INSERTIONS_DF$START[j], INSERTIONS_DF$STOP[j]), IRanges(INSERTIONS_DF$START, INSERTIONS_DF$STOP), type="equal", maxgap = distanceThreshold) ), "ID"]
        INSERTIONS_DF$Samples[j] <- as.character(paste(DUPLICATES, collapse = ";"))
      }
      rm(DUPLICATES)
      
      #B collapse overlapping events
      used <- character()
      for (j in 1:nrow(INSERTIONS_DF)){
        if(!INSERTIONS_DF$ID[j] %in% used){
          # get matches
          matches <- unlist(str_split(INSERTIONS_DF$Samples[j],";"))
          
          if(length(matches)!=0){
            # df
            temp <- INSERTIONS_DF[which(INSERTIONS_DF$ID %in% matches),]
            # find max support CNV and use it
            idx <- which.max(INSERTIONS_DF$MinSupport)
            temp$maxSupSTART[idx] <- temp$START[idx]
            temp$maxSupSTOP[idx] <- temp$STOP[idx]
            temp$MaxSupport[idx] <- max(temp$MinSupport)
            temp$MinSupport[idx] <- min(temp$MinSupport)
            temp$START[idx] <- min(temp$START)
            temp$STOP[idx] <- max(temp$STOP)
            #delete from pool
            used <- append(used, temp$ID)
            temp <- temp[idx,]
            
            #append new
            INSERTIONS_DF_tmp <- rbind(INSERTIONS_DF_tmp,temp)
            
          }else(
            #append current
            INSERTIONS_DF_tmp <- rbind(INSERTIONS_DF_tmp, INSERTIONS_DF[j,])
          )
          
        }
      }
      
      #replace and get new IDs
      INSERTIONS_DF <- INSERTIONS_DF_tmp
      INSERTIONS_DF$ID <- c(1:nrow(INSERTIONS_DF))
      INSERTIONS_DF$ID <- paste(contigNum,"INS",INSERTIONS_DF$ID,sep="_")
      INSERTIONS_DF$Samples <- ""
      
      
      #C BACKTRACKING - FIND IN WHICH SAMPLES THEY ARE
      TMPx <- TMPallCNVs[TMPallCNVs$SVTYPE_MERGED=="INS",]
      INSERTIONS_DF$NumSamples <- 0
      # INSERTIONS_DF$IDs <- ""
      for (j in 1:nrow(INSERTIONS_DF)){
        # which samples overlap with the region?
        MATCHES <- TMPx[subjectHits(IRanges::findOverlaps(IRanges(INSERTIONS_DF$START[j], INSERTIONS_DF$STOP[j]), IRanges(TMPx$START, TMPx$STOP),type="equal",maxgap = distanceThreshold) ),]
        INSERTIONS_DF$Samples[j] <- paste(as.character(unique(MATCHES$SAMPLE)), collapse = ";")
        INSERTIONS_DF$NumSamples[j] <- length(unique(MATCHES$SAMPLE))
        # INSERTIONS_DF$IDs[j] <- paste(as.character(unique(MATCHES$ID)), collapse = ";")
      }
      INSERTIONS_DF <- INSERTIONS_DF[,c("ID","CHROMOSOME","START","STOP","LENGTH_MERGED","SVTYPE_MERGED", "NumSamples", "Samples")]
      INSERTIONS_DF <- INSERTIONS_DF[INSERTIONS_DF$NumSamples>=minSamples,]
    }
    
    
    
    
    
    
    #######################################################################################################################
    ##INVERSIONS
    INVERSIONS_DF <- data.frame(ID=character(),CHROMOSOME=character(),START=integer(), STOP=integer(),
                                LENGTH_MERGED=integer(),SVTYPE_MERGED=character(), MinSupport=integer(),MaxSupport=integer(),
                                maxSupSTART=integer(),maxSupSTOP=integer(),
                                Samples=character()
    )
    INVERSIONS_DF_tmp <- INVERSIONS_DF
    
    
    if(as.integer(max(INVERSIONS))>0){
      
      #A get regions boundaries
      for(level in seq(as.integer(max(INVERSIONS)),threshold,by=-1) ){
        # boundaries https://stackoverflow.com/questions/29184297/finding-the-start-and-stop-indices-in-sequence-in-r
        y <- which(INVERSIONS>=level)
        START <- y[!(y-1) %in% y]
        STOP <- y[!(y+1) %in% y]
        # dataframe filling
        ID <- c(1:length(STOP)) 
        DF_temp <- as.data.frame(cbind(ID,START, STOP))
        DF_temp$LENGTH_MERGED <- DF_temp$STOP - DF_temp$START+1
        DF_temp$MinSupport <- level
        DF_temp$MaxSupport <- 0
        DF_temp$maxSupSTART <- 0
        DF_temp$maxSupSTOP <- 0
        DF_temp$Samples <- ""
        # append 
        INVERSIONS_DF <- rbind(INVERSIONS_DF,DF_temp)
      }
      
      # modify DF columns
      INVERSIONS_DF$CHROMOSOME <- refHeader
      INVERSIONS_DF$ID <- c(1:nrow(INVERSIONS_DF))
      INVERSIONS_DF$ID <- paste(contigNum,"INV",INVERSIONS_DF$ID,sep="_")
      INVERSIONS_DF$SVTYPE_MERGED <- "INV"
      
    }
    
    #B find overlapping events
    if(nrow(INVERSIONS_DF)>0){
      
      for (j in 1:nrow(INVERSIONS_DF)){
        # https://stackoverflow.com/questions/37754509/finding-overlapping-intervals
        DUPLICATES <- INVERSIONS_DF[subjectHits(IRanges::findOverlaps(IRanges(INVERSIONS_DF$START[j], INVERSIONS_DF$STOP[j]), IRanges(INVERSIONS_DF$START, INVERSIONS_DF$STOP), type="equal", maxgap = distanceThreshold) ), "ID"]
        INVERSIONS_DF$Samples[j] <- as.character(paste(DUPLICATES, collapse = ";"))
      }
      rm(DUPLICATES)
      
      #B collapse overlapping events
      used <- character()
      for (j in 1:nrow(INVERSIONS_DF)){
        if(!INVERSIONS_DF$ID[j] %in% used){
          # get matches
          matches <- unlist(str_split(INVERSIONS_DF$Samples[j],";"))
          
          if(length(matches)!=0){
            # df
            temp <- INVERSIONS_DF[which(INVERSIONS_DF$ID %in% matches),]
            # find max support CNV and use it
            idx <- which.max(INVERSIONS_DF$MinSupport)
            temp$maxSupSTART[idx] <- temp$START[idx]
            temp$maxSupSTOP[idx] <- temp$STOP[idx]
            temp$MaxSupport[idx] <- max(temp$MinSupport)
            temp$MinSupport[idx] <- min(temp$MinSupport)
            temp$START[idx] <- min(temp$START)
            temp$STOP[idx] <- max(temp$STOP)
            #delete from pool
            used <- append(used, temp$ID)
            temp <- temp[idx,]
            
            #append new
            INVERSIONS_DF_tmp <- rbind(INVERSIONS_DF_tmp,temp)
            
          }else(
            #append current
            INVERSIONS_DF_tmp <- rbind(INVERSIONS_DF_tmp, INVERSIONS_DF[j,])
          )
          
        }
      }
      
      #replace and get new IDs
      INVERSIONS_DF <- INVERSIONS_DF_tmp
      INVERSIONS_DF$ID <- c(1:nrow(INVERSIONS_DF))
      INVERSIONS_DF$ID <- paste(contigNum,"INV",INVERSIONS_DF$ID,sep="_")
      INVERSIONS_DF$Samples <- ""
      
      
      #C BACKTRACKING - FIND IN WHICH SAMPLES THEY ARE
      TMPx <- TMPallCNVs[TMPallCNVs$SVTYPE_MERGED=="INV",]
      INVERSIONS_DF$NumSamples <- 0
      # INVERSIONS_DF$IDs <- ""
      for (j in 1:nrow(INVERSIONS_DF)){
        # which samples overlap with the region?
        MATCHES <- TMPx[subjectHits(IRanges::findOverlaps(IRanges(INVERSIONS_DF$START[j], INVERSIONS_DF$STOP[j]), IRanges(TMPx$START, TMPx$STOP),type="equal",maxgap = distanceThreshold) ),]
        INVERSIONS_DF$Samples[j] <- paste(as.character(unique(MATCHES$SAMPLE)), collapse = ";")
        INVERSIONS_DF$NumSamples[j] <- length(unique(MATCHES$SAMPLE))
        # INVERSIONS_DF$IDs[j] <- paste(as.character(unique(MATCHES$ID)), collapse = ";")
      }
      INVERSIONS_DF <- INVERSIONS_DF[,c("ID","CHROMOSOME","START","STOP","LENGTH_MERGED","SVTYPE_MERGED", "NumSamples", "Samples")]
      INVERSIONS_DF <- INVERSIONS_DF[INVERSIONS_DF$NumSamples>=minSamples,]
    }
    
    
    #######################################################################################################################
    # GET TOGETHER # WHAT ABOUT LARGE FILES?
    SHORT_RESULTS <- rbindlist(list(SHORT_RESULTS, DELETIONS_DF, DUPLICATIONS_DF, INSERTIONS_DF, INVERSIONS_DF), fill=TRUE)
    
  }
  
  #######################################################################################################################
  #sort
  SHORT_RESULTS <- SHORT_RESULTS[order(SHORT_RESULTS$CHROMOSOME, SHORT_RESULTS$START),]
  SHORT_RESULTS <- SHORT_RESULTS[,c("ID","CHROMOSOME","START","STOP","LENGTH_MERGED","SVTYPE_MERGED", "NumSamples", "Samples")]
  
  #write TABLE
  dir.create(file.path(getwd(),dirname(output)), recursive = TRUE,showWarnings = FALSE)
  write.table(SHORT_RESULTS, file=output, sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE )
  
  
  
}


# develop and test
# args <- c("results/merged_procaryaSV/coverage20.procaryaSV_callers_merge.tsv","results/merged_procaryaSV/coverage20.procaryaSV_venn.png","results/merged_procaryaSV/coverage20.procaryaSV_sv_types.png","results/references/FN433596.fasta","10","NA","2","2000","results/lumpy/coverage20/coverage20.vcf","results/delly2/coverage20/coverage20.vcf","results/cnvnator/coverage20/coverage20.vcf","results/pindel/coverage20/coverage20.vcf","results/cnproscan/coverage20/coverage20.vcf")

# setwd("/home/rj/4TB/PHD_TESTING/artificial_ProcaryaSV/")

#run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)