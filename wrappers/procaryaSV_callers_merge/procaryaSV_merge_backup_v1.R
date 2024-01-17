library(data.table)
library(stringr)
library(seqinr)

## TO DO:
# nastavitelny prah pro min a max delku SV, napr delly a lumpy maji nektere FP velmi dlouhe, max_sv_length - pozor definovane pozdeji, smazat!
# pzor INSERCE maji nulovou delku - nesmazat je
# pridat kontrolu podminky zda nejaky vstupni VCF chybi, co pak?
# parametr min. callers?
# overlap parametr?
# je lepsi nepredpokladat stale stejny header - nahradit header-FALSE?, ale zase pak kontrolovat co je v jakem sloupci, takze nechat header TRUE
# anotace genbank filem - udelat pres HOMER annotatePeaks
## POZOR foverlaps reportuje vsechny i dilci overlapy, proto Subevents~=Pocet overlapu s danym callerem
# ORIGRECORDS <- ORIGRECORDS[-which(ORIGRECORDS$SVTYPE_MERGED!=ORIGRECORDS$SVTYPE),] # popremyslet
# can LUMPY detects INS and INV a co dalsi?????
# pridat VENNUV DIAGRAM

# kontrola SUBTYPE TYPE a ne jen souradnice, nekdy je uvnitr DELECE mnoho INV apod, co s tim

# CO SE MUZE STAT:
# uplne chybi VCF
# VCF nechybi ale je prazdne

run_all <- function(args){
  # arguments
  output <- args[1]
  output_short <- args[2]
  fastaFile <- args[3] 
  min_sv_length <- args[4] #default 1 
  max_sv_length <- args[5] #default refLength/5
  vcf_files <- args[6:length(args)]
  
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
  
  
  ## PROCESS CHROMOSOMES
  numberOfFastaContigs
  
  ## for single contig/chromosome
  if(numberOfFastaContigs == 1){
    refHeader <- fastaContiglist[1]
    refLength <- referenceLengthList[1]
    
    if(is.na(min_sv_length)){min_sv_length <- 1} #default
    if(is.na(max_sv_length)){max_sv_length <- refLength/5} #default
    
    # DELLY2
    if(!identical(delly, character(0))){
      TMPdellyDF <- dellyDF[dellyDF$`#CHROM`==refHeader,]
      TMPdellyDF <- TMPdellyDF[TMPdellyDF$LEN<=max_sv_length,] #delete SV calls which are above 1/5 of the genome length, which is likely False Positive
      TMPdellyDF <- TMPdellyDF[TMPdellyDF$LEN>=min_sv_length,]
      dellyVEC_DUP <- rep(0,refLength)
      dellyVEC_DEL <- rep(0,refLength)
      dellyVEC_INS <- rep(0,refLength)
      dellyVEC_INV <- rep(0,refLength)
      for (j in 1:nrow(TMPdellyDF)){
        if(TMPdellyDF$SVTYPE[j]=="DUP"){dellyVEC_DUP[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
        if(TMPdellyDF$SVTYPE[j]=="DEL"){dellyVEC_DEL[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
        if(TMPdellyDF$SVTYPE[j]=="INS"){dellyVEC_INS[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
        if(TMPdellyDF$SVTYPE[j]=="INV"){dellyVEC_INV[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
      }
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
      TMPlumpyDF <- TMPlumpyDF[TMPlumpyDF$LEN<=max_sv_length,] #delete SV calls which are above 1/5 of the genome length, which is likely False Positive
      TMPlumpyDF <- TMPlumpyDF[TMPlumpyDF$LEN>=min_sv_length,]
      lumpyVEC_DUP <- rep(0,refLength)
      lumpyVEC_DEL <- rep(0,refLength)
      for (j in 1:nrow(TMPlumpyDF)){
        if(TMPlumpyDF$SVTYPE[j]=="DUP"){lumpyVEC_DUP[TMPlumpyDF$POS[j]:TMPlumpyDF$END[j]] <- 1}
        if(TMPlumpyDF$SVTYPE[j]=="DEL"){lumpyVEC_DEL[TMPlumpyDF$POS[j]:TMPlumpyDF$END[j]] <- 1}
      }
    }else{
      lumpyVEC_DUP <- rep(0,refLength)
      lumpyVEC_DEL <- rep(0,refLength)
      TMPlumpyDF <- lumpyDF
    }
    
    # cnvnator
    if(!identical(cnvnator, character(0))){
      TMPcnvnatorDF <- cnvnatorDF[cnvnatorDF$`#CHROM`==refHeader,]
      TMPcnvnatorDF <- TMPcnvnatorDF[TMPcnvnatorDF$LEN<=max_sv_length,] #delete SV calls which are above 1/5 of the genome length, which is likely False Positive
      TMPcnvnatorDF <- TMPcnvnatorDF[TMPcnvnatorDF$LEN>=min_sv_length,]
      cnvnatorVEC_DUP <- rep(0,refLength)
      cnvnatorVEC_DEL <- rep(0,refLength)
      for (j in 1:nrow(TMPcnvnatorDF)){
        if(TMPcnvnatorDF$SVTYPE[j]=="DUP"){cnvnatorVEC_DUP[TMPcnvnatorDF$POS[j]:TMPcnvnatorDF$END[j]] <- 1}
        if(TMPcnvnatorDF$SVTYPE[j]=="DEL"){cnvnatorVEC_DEL[TMPcnvnatorDF$POS[j]:TMPcnvnatorDF$END[j]] <- 1}
      }
    }else{
      cnvnatorVEC_DUP <- rep(0,refLength)
      cnvnatorVEC_DEL <- rep(0,refLength)
      TMPcnvnatorDF <- cnvnatorDF
    }
    
    # pindel
    if(!identical(pindel, character(0))){
      TMPpindelDF <- pindelDF[pindelDF$`#CHROM`==refHeader,]
      TMPpindelDF <- TMPpindelDF[TMPpindelDF$LEN<=max_sv_length,] #delete SV calls which are above 1/5 of the genome length, which is likely False Positive
      TMPpindelDF <- TMPpindelDF[TMPpindelDF$LEN>=min_sv_length,]
      pindelVEC_DUP <- rep(0,refLength)
      pindelVEC_DEL <- rep(0,refLength)
      pindelVEC_INS <- rep(0,refLength)
      pindelVEC_INV <- rep(0,refLength)
      # ingoring RPL - replacement https://www.biostars.org/p/113765/
      for (j in 1:nrow(TMPpindelDF)){
        if(TMPpindelDF$SVTYPE[j]=="DUP"){pindelVEC_DUP[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
        if(TMPpindelDF$SVTYPE[j]=="DEL"){pindelVEC_DEL[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
        if(TMPpindelDF$SVTYPE[j]=="INS"){pindelVEC_INS[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
        if(TMPpindelDF$SVTYPE[j]=="INV"){pindelVEC_INV[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
      }
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
      TMPcnproscanDF <- TMPcnproscanDF[TMPcnproscanDF$LEN<=max_sv_length,] #delete SV calls which are above 1/5 of the genome length, which is likely False Positive
      TMPcnproscanDF <- TMPcnproscanDF[TMPcnproscanDF$LEN>=min_sv_length,]
      cnproscanVEC_DUP <- rep(0,refLength)
      cnproscanVEC_DEL <- rep(0,refLength)
      # ingoring RPL - replacement https://www.biostars.org/p/113765/
      for (j in 1:nrow(TMPcnproscanDF)){
        if(TMPcnproscanDF$SVTYPE[j]=="DUP"){cnproscanVEC_DUP[TMPcnproscanDF$POS[j]:TMPcnproscanDF$END[j]] <- 1}
        if(TMPcnproscanDF$SVTYPE[j]=="DEL"){cnproscanVEC_DEL[TMPcnproscanDF$POS[j]:TMPcnproscanDF$END[j]] <- 1}
      }
    }else{
      cnproscanVEC_DUP <- rep(0,refLength)
      cnproscanVEC_DEL <- rep(0,refLength)
      TMPcnproscanDF <- cnproscanDF
    }
    
    ## data.table preparation
    data.table::setDT(TMPdellyDF)
    data.table::setDT(TMPlumpyDF)
    data.table::setDT(TMPcnvnatorDF)
    data.table::setDT(TMPpindelDF)
    data.table::setDT(TMPcnproscanDF)
    data.table::setkey(TMPdellyDF,"#CHROM","POS","END")
    data.table::setkey(TMPlumpyDF,"#CHROM","POS","END")
    data.table::setkey(TMPcnvnatorDF,"#CHROM","POS","END")
    data.table::setkey(TMPpindelDF,"#CHROM","POS","END")
    data.table::setkey(TMPcnproscanDF,"#CHROM","POS","END")
    
    ## CREATE VECTORS
    DELETIONS <- dellyVEC_DEL+lumpyVEC_DEL+cnvnatorVEC_DEL+pindelVEC_DEL+cnproscanVEC_DEL
    DUPLICATIONS <- dellyVEC_DUP+lumpyVEC_DUP+cnvnatorVEC_DUP+pindelVEC_DUP+cnproscanVEC_DUP
    INSERTIONS <- pindelVEC_INS+dellyVEC_INS
    INVERSIONS <- pindelVEC_INV+dellyVEC_INV
    
    #######################################################################################################################
    ## DELETIONS
    #GET START and STOP of EVENTS - https://stackoverflow.com/questions/29184297/finding-the-start-and-stop-indices-in-sequence-in-r
    y <- which(DELETIONS>0)
    START <- y[!(y-1) %in% y]
    STOP <- y[!(y+1) %in% y]
    DELETIONS_DF <- as.data.frame(cbind(START, STOP))
    DELETIONS_DF$LEN <- DELETIONS_DF$STOP - DELETIONS_DF$START+1
    DELETIONS_DF$NumOfCallers <- 0
    DELETIONS_DF <- DELETIONS_DF[DELETIONS_DF$LEN>min_sv_length,]
    DELETIONS_DF$Callers <- ""
    
    DELETIONS_DF$cnproscanPortion <- 0
    DELETIONS_DF$delly2Portion <- 0
    DELETIONS_DF$lumpyPortion <- 0
    DELETIONS_DF$cnvnatorPortion <- 0
    DELETIONS_DF$pindelPortion <- 0
    
    DELETIONS_DF$cnproscanSubevents <- 0
    DELETIONS_DF$delly2Subevents <- 0
    DELETIONS_DF$lumpySubevents <- 0
    DELETIONS_DF$cnvnatorSubevents <- 0
    DELETIONS_DF$pindelSubevents <- 0
    
    # GET INFO
    for (j in 1:nrow(DELETIONS_DF)){
      DELETIONS_DF$NumOfCallers[j] <- max(DELETIONS[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])])
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
    DELETIONS_DF <- DELETIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                    "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                    "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )]
    colnames(DELETIONS_DF) <- c("CHROMOSOME","START","STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )
    
    ## get back original records foverlaps
    data.table::setDT(DELETIONS_DF)
    data.table::setkey(DELETIONS_DF,"CHROMOSOME","START","STOP")
    
    ORIGRECORDS_delly <- data.table::foverlaps(TMPdellyDF, DELETIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    ORIGRECORDS_lumpy <- data.table::foverlaps(TMPlumpyDF, DELETIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    ORIGRECORDS_cnvnator <- data.table::foverlaps(TMPcnvnatorDF, DELETIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    ORIGRECORDS_pindel <- data.table::foverlaps(TMPpindelDF, DELETIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    ORIGRECORDS_cnproscan <- data.table::foverlaps(TMPcnproscanDF, DELETIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    
    ORIGRECORDS_DELETIONS <- rbind(ORIGRECORDS_delly,ORIGRECORDS_lumpy,ORIGRECORDS_cnvnator,ORIGRECORDS_pindel,ORIGRECORDS_cnproscan,fill=TRUE)
    
    
    colnames(ORIGRECORDS_DELETIONS)[colnames(ORIGRECORDS_DELETIONS) == "#CHROM"] <- "CHROMOSOME"
    colnames(ORIGRECORDS_DELETIONS)[colnames(ORIGRECORDS_DELETIONS) == "POS"] <- "START_EVENT"
    colnames(ORIGRECORDS_DELETIONS)[colnames(ORIGRECORDS_DELETIONS) == "END"] <- "END_EVENT"
    colnames(ORIGRECORDS_DELETIONS)[colnames(ORIGRECORDS_DELETIONS) == "CALLER"] <- "EVENT_CALLER"
    colnames(ORIGRECORDS_DELETIONS)[colnames(ORIGRECORDS_DELETIONS) == "LEN"] <- "EVENT_LENGTH"
    ORIGRECORDS_DELETIONS$SEPARATOR <- ""
    ORIGRECORDS_DELETIONS <- ORIGRECORDS_DELETIONS[,c("CHROMOSOME","START" ,"STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                                      "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                                      "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents",
                                                      "SEPARATOR",
                                                      "EVENT_CALLER", "START_EVENT","END_EVENT","SVTYPE",              
                                                      "VCF_INFO","EVENT_LENGTH" )]
    
    ORIGRECORDS_DELETIONS <- unique(ORIGRECORDS_DELETIONS)
    
    #######################################################################################################################
    ## DUPLICATIONS
    y <- which(DUPLICATIONS>0)
    START <- y[!(y-1) %in% y]
    STOP <- y[!(y+1) %in% y]
    DUPLICATIONS_DF <- as.data.frame(cbind(START, STOP))
    DUPLICATIONS_DF$LEN <- DUPLICATIONS_DF$STOP - DUPLICATIONS_DF$START+1
    DUPLICATIONS_DF$NumOfCallers <- 0
    DUPLICATIONS_DF <- DUPLICATIONS_DF[DUPLICATIONS_DF$LEN>min_sv_length,]
    DUPLICATIONS_DF$Callers <- ""
    
    DUPLICATIONS_DF$cnproscanPortion <- 0
    DUPLICATIONS_DF$delly2Portion <- 0
    DUPLICATIONS_DF$lumpyPortion <- 0
    DUPLICATIONS_DF$cnvnatorPortion <- 0
    DUPLICATIONS_DF$pindelPortion <- 0
    
    DUPLICATIONS_DF$cnproscanSubevents <- 0
    DUPLICATIONS_DF$delly2Subevents <- 0
    DUPLICATIONS_DF$lumpySubevents <- 0
    DUPLICATIONS_DF$cnvnatorSubevents <- 0
    DUPLICATIONS_DF$pindelSubevents <- 0
    
    # GET INFO
    for (j in 1:nrow(DUPLICATIONS_DF)){
      DUPLICATIONS_DF$NumOfCallers[j] <- max(DUPLICATIONS[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])])
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
    DUPLICATIONS_DF <- DUPLICATIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                          "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                          "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )]
    colnames(DUPLICATIONS_DF) <- c("CHROMOSOME","START","STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                   "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                   "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents"  )
    
    ## get back original records foverlaps
    data.table::setDT(DUPLICATIONS_DF)
    data.table::setkey(DUPLICATIONS_DF,"CHROMOSOME","START","STOP")
    
    ORIGRECORDS_delly <- data.table::foverlaps(TMPdellyDF, DUPLICATIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    ORIGRECORDS_lumpy <- data.table::foverlaps(TMPlumpyDF, DUPLICATIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    ORIGRECORDS_cnvnator <- data.table::foverlaps(TMPcnvnatorDF, DUPLICATIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    ORIGRECORDS_pindel <- data.table::foverlaps(TMPpindelDF, DUPLICATIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    ORIGRECORDS_cnproscan <- data.table::foverlaps(TMPcnproscanDF, DUPLICATIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    
    ORIGRECORDS_DUPLICATIONS <- rbind(ORIGRECORDS_delly,ORIGRECORDS_lumpy,ORIGRECORDS_cnvnator,ORIGRECORDS_pindel,ORIGRECORDS_cnproscan,fill=TRUE)
    
    
    colnames(ORIGRECORDS_DUPLICATIONS)[colnames(ORIGRECORDS_DUPLICATIONS) == "#CHROM"] <- "CHROMOSOME"
    colnames(ORIGRECORDS_DUPLICATIONS)[colnames(ORIGRECORDS_DUPLICATIONS) == "POS"] <- "START_EVENT"
    colnames(ORIGRECORDS_DUPLICATIONS)[colnames(ORIGRECORDS_DUPLICATIONS) == "END"] <- "END_EVENT"
    colnames(ORIGRECORDS_DUPLICATIONS)[colnames(ORIGRECORDS_DUPLICATIONS) == "CALLER"] <- "EVENT_CALLER"
    colnames(ORIGRECORDS_DUPLICATIONS)[colnames(ORIGRECORDS_DUPLICATIONS) == "LEN"] <- "EVENT_LENGTH"
    ORIGRECORDS_DUPLICATIONS$SEPARATOR <- ""
    ORIGRECORDS_DUPLICATIONS <- ORIGRECORDS_DUPLICATIONS[,c("CHROMOSOME","START" ,"STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                                            "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                                            "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents",
                                                            "SEPARATOR", "EVENT_CALLER", "START_EVENT","END_EVENT","SVTYPE",              
                                                            "VCF_INFO","EVENT_LENGTH" )]
    
    ORIGRECORDS_DUPLICATIONS <- unique(ORIGRECORDS_DUPLICATIONS)
    
    #######################################################################################################################
    ## INSERTIONS
    y <- which(INSERTIONS>0)
    START <- y[!(y-1) %in% y]
    STOP <- y[!(y+1) %in% y]
    INSERTIONS_DF <- as.data.frame(cbind(START, STOP))
    INSERTIONS_DF$LEN <- INSERTIONS_DF$STOP - INSERTIONS_DF$START+1
    INSERTIONS_DF$NumOfCallers <- 0
    # INSERTIONS_DF <- INSERTIONS_DF[INSERTIONS_DF$LEN>min_sv_length,]
    INSERTIONS_DF$Callers <- ""
    
    INSERTIONS_DF$cnproscanPortion <- 0
    INSERTIONS_DF$delly2Portion <- 0
    INSERTIONS_DF$lumpyPortion <- 0
    INSERTIONS_DF$cnvnatorPortion <- 0
    INSERTIONS_DF$pindelPortion <- 0
    
    INSERTIONS_DF$cnproscanSubevents <- 0
    INSERTIONS_DF$delly2Subevents <- 0
    INSERTIONS_DF$lumpySubevents <- 0
    INSERTIONS_DF$cnvnatorSubevents <- 0
    INSERTIONS_DF$pindelSubevents <- 0
    
    # GET INFO
    for (j in 1:nrow(INSERTIONS_DF)){
      INSERTIONS_DF$NumOfCallers[j] <- max(INSERTIONS[(INSERTIONS_DF$START[j]):(INSERTIONS_DF$STOP[j])])
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
    INSERTIONS_DF <- INSERTIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                      "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                      "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )]
    colnames(INSERTIONS_DF) <- c("CHROMOSOME","START","STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                 "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                 "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )
    
    ## get back original records foverlaps
    data.table::setDT(INSERTIONS_DF)
    data.table::setkey(INSERTIONS_DF,"CHROMOSOME","START","STOP")
    
    ORIGRECORDS_delly <- data.table::foverlaps(TMPdellyDF, INSERTIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    ORIGRECORDS_pindel <- data.table::foverlaps(TMPpindelDF, INSERTIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    
    ORIGRECORDS_INSERTIONS <- rbind(ORIGRECORDS_delly,ORIGRECORDS_pindel,fill=TRUE) #ORIGRECORDS_lumpy
    
    colnames(ORIGRECORDS_INSERTIONS)[colnames(ORIGRECORDS_INSERTIONS) == "#CHROM"] <- "CHROMOSOME"
    colnames(ORIGRECORDS_INSERTIONS)[colnames(ORIGRECORDS_INSERTIONS) == "POS"] <- "START_EVENT"
    colnames(ORIGRECORDS_INSERTIONS)[colnames(ORIGRECORDS_INSERTIONS) == "END"] <- "END_EVENT"
    colnames(ORIGRECORDS_INSERTIONS)[colnames(ORIGRECORDS_INSERTIONS) == "CALLER"] <- "EVENT_CALLER"
    colnames(ORIGRECORDS_INSERTIONS)[colnames(ORIGRECORDS_INSERTIONS) == "LEN"] <- "EVENT_LENGTH"
    ORIGRECORDS_INSERTIONS$SEPARATOR <- ""
    ORIGRECORDS_INSERTIONS <- ORIGRECORDS_INSERTIONS[,c("CHROMOSOME","START" ,"STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                                        "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                                        "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents",
                                                        "SEPARATOR",
                                                        "EVENT_CALLER", "START_EVENT","END_EVENT","SVTYPE",              
                                                        "VCF_INFO","EVENT_LENGTH" )]
    
    ORIGRECORDS_INSERTIONS <- unique(ORIGRECORDS_INSERTIONS)
    
    #######################################################################################################################
    ##INVERSIONS
    y <- which(INVERSIONS>0)
    START <- y[!(y-1) %in% y]
    STOP <- y[!(y+1) %in% y]
    INVERSIONS_DF <- as.data.frame(cbind(START, STOP))
    INVERSIONS_DF$LEN <- INVERSIONS_DF$STOP - INVERSIONS_DF$START+1
    INVERSIONS_DF$NumOfCallers <- 0
    INVERSIONS_DF <- INVERSIONS_DF[INVERSIONS_DF$LEN>min_sv_length,]
    INVERSIONS_DF$Callers <- ""
    
    INVERSIONS_DF$cnproscanPortion <- 0
    INVERSIONS_DF$delly2Portion <- 0
    INVERSIONS_DF$lumpyPortion <- 0
    INVERSIONS_DF$cnvnatorPortion <- 0
    INVERSIONS_DF$pindelPortion <- 0
    
    INVERSIONS_DF$cnproscanSubevents <- 0
    INVERSIONS_DF$delly2Subevents <- 0
    INVERSIONS_DF$lumpySubevents <- 0
    INVERSIONS_DF$cnvnatorSubevents <- 0
    INVERSIONS_DF$pindelSubevents <- 0
    
    # GET INFO
    for (j in 1:nrow(INVERSIONS_DF)){
      INVERSIONS_DF$NumOfCallers[j] <- max(INVERSIONS[(INVERSIONS_DF$START[j]):(INVERSIONS_DF$STOP[j])])
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
    INVERSIONS_DF <- INVERSIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                      "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                      "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents")]
    colnames(INVERSIONS_DF) <- c("CHROMOSOME","START","STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                 "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                 "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents")
    
    ## get back original records foverlaps
    data.table::setDT(INVERSIONS_DF)
    data.table::setkey(INVERSIONS_DF,"CHROMOSOME","START","STOP")
    
    ORIGRECORDS_delly <- data.table::foverlaps(TMPdellyDF, INVERSIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    ORIGRECORDS_pindel <- data.table::foverlaps(TMPpindelDF, INVERSIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
    
    ORIGRECORDS_INVERSIONS <- rbind(ORIGRECORDS_delly,ORIGRECORDS_pindel,fill=TRUE) #ORIGRECORDS_lumpy
    
    colnames(ORIGRECORDS_INVERSIONS)[colnames(ORIGRECORDS_INVERSIONS) == "#CHROM"] <- "CHROMOSOME"
    colnames(ORIGRECORDS_INVERSIONS)[colnames(ORIGRECORDS_INVERSIONS) == "POS"] <- "START_EVENT"
    colnames(ORIGRECORDS_INVERSIONS)[colnames(ORIGRECORDS_INVERSIONS) == "END"] <- "END_EVENT"
    colnames(ORIGRECORDS_INVERSIONS)[colnames(ORIGRECORDS_INVERSIONS) == "CALLER"] <- "EVENT_CALLER"
    colnames(ORIGRECORDS_INVERSIONS)[colnames(ORIGRECORDS_INVERSIONS) == "LEN"] <- "EVENT_LENGTH"
    ORIGRECORDS_INVERSIONS$SEPARATOR <- ""
    ORIGRECORDS_INVERSIONS <- ORIGRECORDS_INVERSIONS[,c("CHROMOSOME","START" ,"STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                                        "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                                        "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents",
                                                        "SEPARATOR", "EVENT_CALLER", "START_EVENT","END_EVENT","SVTYPE",              
                                                        "VCF_INFO","EVENT_LENGTH" )]
    
    ORIGRECORDS_INVERSIONS <- unique(ORIGRECORDS_INVERSIONS)
    
    # GET TOGETHER # WHAT ABOUT LARGE FILES?
    RESULTS <- rbindlist(list(ORIGRECORDS_DELETIONS,ORIGRECORDS_DUPLICATIONS,ORIGRECORDS_INSERTIONS,ORIGRECORDS_INVERSIONS),fill=TRUE)
    SHORT_RESULTS <- rbindlist(list(DELETIONS_DF, DUPLICATIONS_DF, INSERTIONS_DF, INVERSIONS_DF),fill=TRUE)
  }
  
  
  # for more contigs/chromosomes
  if(numberOfFastaContigs  > 1){
    
    RESULTS <- data.table()
    SHORT_RESULTS <- data.table()
    
    for (i in 1:numberOfFastaContigs){
      
      refHeader <- fastaContiglist[i]
      refLength <- referenceLengthList[i]
      
      if(is.na(min_sv_length)){min_sv_length <- 1} #default
      if(is.na(max_sv_length)){max_sv_length <- refLength/5} #default
      
      # DELLY2
      if(!identical(delly, character(0))){
        TMPdellyDF <- dellyDF[dellyDF$`#CHROM`==refHeader,]
        TMPdellyDF <- TMPdellyDF[TMPdellyDF$LEN<=max_sv_length,] #delete SV calls which are above 1/5 of the genome length, which is likely False Positive
        TMPdellyDF <- TMPdellyDF[TMPdellyDF$LEN>=min_sv_length,]
        dellyVEC_DUP <- rep(0,refLength)
        dellyVEC_DEL <- rep(0,refLength)
        dellyVEC_INS <- rep(0,refLength)
        dellyVEC_INV <- rep(0,refLength)
        for (j in 1:nrow(TMPdellyDF)){
          if(TMPdellyDF$SVTYPE[j]=="DUP"){dellyVEC_DUP[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
          if(TMPdellyDF$SVTYPE[j]=="DEL"){dellyVEC_DEL[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
          if(TMPdellyDF$SVTYPE[j]=="INS"){dellyVEC_INS[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
          if(TMPdellyDF$SVTYPE[j]=="INV"){dellyVEC_INV[TMPdellyDF$POS[j]:TMPdellyDF$END[j]] <- 1}
        }
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
        TMPlumpyDF <- TMPlumpyDF[TMPlumpyDF$LEN<=max_sv_length,] #delete SV calls which are above 1/5 of the genome length, which is likely False Positive
        TMPlumpyDF <- TMPlumpyDF[TMPlumpyDF$LEN>=min_sv_length,]
        lumpyVEC_DUP <- rep(0,refLength)
        lumpyVEC_DEL <- rep(0,refLength)
        for (j in 1:nrow(TMPlumpyDF)){
          if(TMPlumpyDF$SVTYPE[j]=="DUP"){lumpyVEC_DUP[TMPlumpyDF$POS[j]:TMPlumpyDF$END[j]] <- 1}
          if(TMPlumpyDF$SVTYPE[j]=="DEL"){lumpyVEC_DEL[TMPlumpyDF$POS[j]:TMPlumpyDF$END[j]] <- 1}
        }
      }else{
        lumpyVEC_DUP <- rep(0,refLength)
        lumpyVEC_DEL <- rep(0,refLength)
        TMPlumpyDF <- lumpyDF
      }
      
      # cnvnator
      if(!identical(cnvnator, character(0))){
        TMPcnvnatorDF <- cnvnatorDF[cnvnatorDF$`#CHROM`==refHeader,]
        TMPcnvnatorDF <- TMPcnvnatorDF[TMPcnvnatorDF$LEN<=max_sv_length,] #delete SV calls which are above 1/5 of the genome length, which is likely False Positive
        TMPcnvnatorDF <- TMPcnvnatorDF[TMPcnvnatorDF$LEN>=min_sv_length,]
        cnvnatorVEC_DUP <- rep(0,refLength)
        cnvnatorVEC_DEL <- rep(0,refLength)
        for (j in 1:nrow(TMPcnvnatorDF)){
          if(TMPcnvnatorDF$SVTYPE[j]=="DUP"){cnvnatorVEC_DUP[TMPcnvnatorDF$POS[j]:TMPcnvnatorDF$END[j]] <- 1}
          if(TMPcnvnatorDF$SVTYPE[j]=="DEL"){cnvnatorVEC_DEL[TMPcnvnatorDF$POS[j]:TMPcnvnatorDF$END[j]] <- 1}
        }
      }else{
        cnvnatorVEC_DUP <- rep(0,refLength)
        cnvnatorVEC_DEL <- rep(0,refLength)
        TMPcnvnatorDF <- cnvnatorDF
      }
      
      # pindel
      if(!identical(pindel, character(0))){
        TMPpindelDF <- pindelDF[pindelDF$`#CHROM`==refHeader,]
        TMPpindelDF <- TMPpindelDF[TMPpindelDF$LEN<=max_sv_length,] #delete SV calls which are above 1/5 of the genome length, which is likely False Positive
        TMPpindelDF <- TMPpindelDF[TMPpindelDF$LEN>=min_sv_length,]
        pindelVEC_DUP <- rep(0,refLength)
        pindelVEC_DEL <- rep(0,refLength)
        pindelVEC_INS <- rep(0,refLength)
        pindelVEC_INV <- rep(0,refLength)
        # ingoring RPL - replacement https://www.biostars.org/p/113765/
        for (j in 1:nrow(TMPpindelDF)){
          if(TMPpindelDF$SVTYPE[j]=="DUP"){pindelVEC_DUP[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
          if(TMPpindelDF$SVTYPE[j]=="DEL"){pindelVEC_DEL[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
          if(TMPpindelDF$SVTYPE[j]=="INS"){pindelVEC_INS[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
          if(TMPpindelDF$SVTYPE[j]=="INV"){pindelVEC_INV[TMPpindelDF$POS[j]:TMPpindelDF$END[j]] <- 1}
        }
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
        TMPcnproscanDF <- TMPcnproscanDF[TMPcnproscanDF$LEN<=max_sv_length,] #delete SV calls which are above 1/5 of the genome length, which is likely False Positive
        TMPcnproscanDF <- TMPcnproscanDF[TMPcnproscanDF$LEN>=min_sv_length,]
        cnproscanVEC_DUP <- rep(0,refLength)
        cnproscanVEC_DEL <- rep(0,refLength)
        # ingoring RPL - replacement https://www.biostars.org/p/113765/
        for (j in 1:nrow(TMPcnproscanDF)){
          if(TMPcnproscanDF$SVTYPE[j]=="DUP"){cnproscanVEC_DUP[TMPcnproscanDF$POS[j]:TMPcnproscanDF$END[j]] <- 1}
          if(TMPcnproscanDF$SVTYPE[j]=="DEL"){cnproscanVEC_DEL[TMPcnproscanDF$POS[j]:TMPcnproscanDF$END[j]] <- 1}
        }
      }else{
        cnproscanVEC_DUP <- rep(0,refLength)
        cnproscanVEC_DEL <- rep(0,refLength)
        TMPcnproscanDF <- cnproscanDF
      }
      
      ## data.table preparation
      data.table::setDT(TMPdellyDF)
      data.table::setDT(TMPlumpyDF)
      data.table::setDT(TMPcnvnatorDF)
      data.table::setDT(TMPpindelDF)
      data.table::setDT(TMPcnproscanDF)
      data.table::setkey(TMPdellyDF,"#CHROM","POS","END")
      data.table::setkey(TMPlumpyDF,"#CHROM","POS","END")
      data.table::setkey(TMPcnvnatorDF,"#CHROM","POS","END")
      data.table::setkey(TMPpindelDF,"#CHROM","POS","END")
      data.table::setkey(TMPcnproscanDF,"#CHROM","POS","END")
      
      ## CREATE VECTORS
      DELETIONS <- dellyVEC_DEL+lumpyVEC_DEL+cnvnatorVEC_DEL+pindelVEC_DEL+cnproscanVEC_DEL
      DUPLICATIONS <- dellyVEC_DUP+lumpyVEC_DUP+cnvnatorVEC_DUP+pindelVEC_DUP+cnproscanVEC_DUP
      INSERTIONS <- pindelVEC_INS+dellyVEC_INS
      INVERSIONS <- pindelVEC_INV+dellyVEC_INV
      
      #######################################################################################################################
      ## DELETIONS
      #GET START and STOP of EVENTS - https://stackoverflow.com/questions/29184297/finding-the-start-and-stop-indices-in-sequence-in-r
      y <- which(DELETIONS>0)
      START <- y[!(y-1) %in% y]
      STOP <- y[!(y+1) %in% y]
      DELETIONS_DF <- as.data.frame(cbind(START, STOP))
      DELETIONS_DF$LEN <- DELETIONS_DF$STOP - DELETIONS_DF$START+1
      DELETIONS_DF$NumOfCallers <- 0
      DELETIONS_DF <- DELETIONS_DF[DELETIONS_DF$LEN>min_sv_length,]
      DELETIONS_DF$Callers <- ""
      
      DELETIONS_DF$cnproscanPortion <- 0
      DELETIONS_DF$delly2Portion <- 0
      DELETIONS_DF$lumpyPortion <- 0
      DELETIONS_DF$cnvnatorPortion <- 0
      DELETIONS_DF$pindelPortion <- 0
      
      DELETIONS_DF$cnproscanSubevents <- 0
      DELETIONS_DF$delly2Subevents <- 0
      DELETIONS_DF$lumpySubevents <- 0
      DELETIONS_DF$cnvnatorSubevents <- 0
      DELETIONS_DF$pindelSubevents <- 0
      
      # GET INFO
      for (j in 1:nrow(DELETIONS_DF)){
        DELETIONS_DF$NumOfCallers[j] <- max(DELETIONS[(DELETIONS_DF$START[j]):(DELETIONS_DF$STOP[j])])
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
      DELETIONS_DF <- DELETIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                      "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                      "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )]
      colnames(DELETIONS_DF) <- c("CHROMOSOME","START","STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                  "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                  "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )
      
      ## get back original records foverlaps
      data.table::setDT(DELETIONS_DF)
      data.table::setkey(DELETIONS_DF,"CHROMOSOME","START","STOP")
      
      ORIGRECORDS_delly <- data.table::foverlaps(TMPdellyDF, DELETIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      ORIGRECORDS_lumpy <- data.table::foverlaps(TMPlumpyDF, DELETIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      ORIGRECORDS_cnvnator <- data.table::foverlaps(TMPcnvnatorDF, DELETIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      ORIGRECORDS_pindel <- data.table::foverlaps(TMPpindelDF, DELETIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      ORIGRECORDS_cnproscan <- data.table::foverlaps(TMPcnproscanDF, DELETIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      
      ORIGRECORDS_DELETIONS <- rbind(ORIGRECORDS_delly,ORIGRECORDS_lumpy,ORIGRECORDS_cnvnator,ORIGRECORDS_pindel,ORIGRECORDS_cnproscan,fill=TRUE)
      
      
      colnames(ORIGRECORDS_DELETIONS)[colnames(ORIGRECORDS_DELETIONS) == "#CHROM"] <- "CHROMOSOME"
      colnames(ORIGRECORDS_DELETIONS)[colnames(ORIGRECORDS_DELETIONS) == "POS"] <- "START_EVENT"
      colnames(ORIGRECORDS_DELETIONS)[colnames(ORIGRECORDS_DELETIONS) == "END"] <- "END_EVENT"
      colnames(ORIGRECORDS_DELETIONS)[colnames(ORIGRECORDS_DELETIONS) == "CALLER"] <- "EVENT_CALLER"
      colnames(ORIGRECORDS_DELETIONS)[colnames(ORIGRECORDS_DELETIONS) == "LEN"] <- "EVENT_LENGTH"
      ORIGRECORDS_DELETIONS$SEPARATOR <- ""
      ORIGRECORDS_DELETIONS <- ORIGRECORDS_DELETIONS[,c("CHROMOSOME","START" ,"STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                                        "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                                        "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents",
                                                        "SEPARATOR",
                                                        "EVENT_CALLER", "START_EVENT","END_EVENT","SVTYPE",              
                                                        "VCF_INFO","EVENT_LENGTH" )]
      
      ORIGRECORDS_DELETIONS <- unique(ORIGRECORDS_DELETIONS)
      
      #######################################################################################################################
      ## DUPLICATIONS
      y <- which(DUPLICATIONS>0)
      START <- y[!(y-1) %in% y]
      STOP <- y[!(y+1) %in% y]
      DUPLICATIONS_DF <- as.data.frame(cbind(START, STOP))
      DUPLICATIONS_DF$LEN <- DUPLICATIONS_DF$STOP - DUPLICATIONS_DF$START+1
      DUPLICATIONS_DF$NumOfCallers <- 0
      DUPLICATIONS_DF <- DUPLICATIONS_DF[DUPLICATIONS_DF$LEN>min_sv_length,]
      DUPLICATIONS_DF$Callers <- ""
      
      DUPLICATIONS_DF$cnproscanPortion <- 0
      DUPLICATIONS_DF$delly2Portion <- 0
      DUPLICATIONS_DF$lumpyPortion <- 0
      DUPLICATIONS_DF$cnvnatorPortion <- 0
      DUPLICATIONS_DF$pindelPortion <- 0
      
      DUPLICATIONS_DF$cnproscanSubevents <- 0
      DUPLICATIONS_DF$delly2Subevents <- 0
      DUPLICATIONS_DF$lumpySubevents <- 0
      DUPLICATIONS_DF$cnvnatorSubevents <- 0
      DUPLICATIONS_DF$pindelSubevents <- 0
      
      # GET INFO
      for (j in 1:nrow(DUPLICATIONS_DF)){
        DUPLICATIONS_DF$NumOfCallers[j] <- max(DUPLICATIONS[(DUPLICATIONS_DF$START[j]):(DUPLICATIONS_DF$STOP[j])])
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
      DUPLICATIONS_DF <- DUPLICATIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                            "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                            "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )]
      colnames(DUPLICATIONS_DF) <- c("CHROMOSOME","START","STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                     "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                     "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents"  )
      
      ## get back original records foverlaps
      data.table::setDT(DUPLICATIONS_DF)
      data.table::setkey(DUPLICATIONS_DF,"CHROMOSOME","START","STOP")
      
      ORIGRECORDS_delly <- data.table::foverlaps(TMPdellyDF, DUPLICATIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      ORIGRECORDS_lumpy <- data.table::foverlaps(TMPlumpyDF, DUPLICATIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      ORIGRECORDS_cnvnator <- data.table::foverlaps(TMPcnvnatorDF, DUPLICATIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      ORIGRECORDS_pindel <- data.table::foverlaps(TMPpindelDF, DUPLICATIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      ORIGRECORDS_cnproscan <- data.table::foverlaps(TMPcnproscanDF, DUPLICATIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      
      ORIGRECORDS_DUPLICATIONS <- rbind(ORIGRECORDS_delly,ORIGRECORDS_lumpy,ORIGRECORDS_cnvnator,ORIGRECORDS_pindel,ORIGRECORDS_cnproscan,fill=TRUE)
      
      
      colnames(ORIGRECORDS_DUPLICATIONS)[colnames(ORIGRECORDS_DUPLICATIONS) == "#CHROM"] <- "CHROMOSOME"
      colnames(ORIGRECORDS_DUPLICATIONS)[colnames(ORIGRECORDS_DUPLICATIONS) == "POS"] <- "START_EVENT"
      colnames(ORIGRECORDS_DUPLICATIONS)[colnames(ORIGRECORDS_DUPLICATIONS) == "END"] <- "END_EVENT"
      colnames(ORIGRECORDS_DUPLICATIONS)[colnames(ORIGRECORDS_DUPLICATIONS) == "CALLER"] <- "EVENT_CALLER"
      colnames(ORIGRECORDS_DUPLICATIONS)[colnames(ORIGRECORDS_DUPLICATIONS) == "LEN"] <- "EVENT_LENGTH"
      ORIGRECORDS_DUPLICATIONS$SEPARATOR <- ""
      ORIGRECORDS_DUPLICATIONS <- ORIGRECORDS_DUPLICATIONS[,c("CHROMOSOME","START" ,"STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                                              "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                                              "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents",
                                                              "SEPARATOR", "EVENT_CALLER", "START_EVENT","END_EVENT","SVTYPE",              
                                                              "VCF_INFO","EVENT_LENGTH" )]
      
      ORIGRECORDS_DUPLICATIONS <- unique(ORIGRECORDS_DUPLICATIONS)
      
      #######################################################################################################################
      ## INSERTIONS
      y <- which(INSERTIONS>0)
      START <- y[!(y-1) %in% y]
      STOP <- y[!(y+1) %in% y]
      INSERTIONS_DF <- as.data.frame(cbind(START, STOP))
      INSERTIONS_DF$LEN <- INSERTIONS_DF$STOP - INSERTIONS_DF$START+1
      INSERTIONS_DF$NumOfCallers <- 0
      # INSERTIONS_DF <- INSERTIONS_DF[INSERTIONS_DF$LEN>min_sv_length,]
      INSERTIONS_DF$Callers <- ""
      
      INSERTIONS_DF$cnproscanPortion <- 0
      INSERTIONS_DF$delly2Portion <- 0
      INSERTIONS_DF$lumpyPortion <- 0
      INSERTIONS_DF$cnvnatorPortion <- 0
      INSERTIONS_DF$pindelPortion <- 0
      
      INSERTIONS_DF$cnproscanSubevents <- 0
      INSERTIONS_DF$delly2Subevents <- 0
      INSERTIONS_DF$lumpySubevents <- 0
      INSERTIONS_DF$cnvnatorSubevents <- 0
      INSERTIONS_DF$pindelSubevents <- 0
      
      # GET INFO
      for (j in 1:nrow(INSERTIONS_DF)){
        INSERTIONS_DF$NumOfCallers[j] <- max(INSERTIONS[(INSERTIONS_DF$START[j]):(INSERTIONS_DF$STOP[j])])
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
      INSERTIONS_DF <- INSERTIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                        "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                        "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )]
      colnames(INSERTIONS_DF) <- c("CHROMOSOME","START","STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                   "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                   "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents" )
      
      ## get back original records foverlaps
      data.table::setDT(INSERTIONS_DF)
      data.table::setkey(INSERTIONS_DF,"CHROMOSOME","START","STOP")
      
      ORIGRECORDS_delly <- data.table::foverlaps(TMPdellyDF, INSERTIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      ORIGRECORDS_pindel <- data.table::foverlaps(TMPpindelDF, INSERTIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      
      ORIGRECORDS_INSERTIONS <- rbind(ORIGRECORDS_delly,ORIGRECORDS_pindel,fill=TRUE) #ORIGRECORDS_lumpy
      
      colnames(ORIGRECORDS_INSERTIONS)[colnames(ORIGRECORDS_INSERTIONS) == "#CHROM"] <- "CHROMOSOME"
      colnames(ORIGRECORDS_INSERTIONS)[colnames(ORIGRECORDS_INSERTIONS) == "POS"] <- "START_EVENT"
      colnames(ORIGRECORDS_INSERTIONS)[colnames(ORIGRECORDS_INSERTIONS) == "END"] <- "END_EVENT"
      colnames(ORIGRECORDS_INSERTIONS)[colnames(ORIGRECORDS_INSERTIONS) == "CALLER"] <- "EVENT_CALLER"
      colnames(ORIGRECORDS_INSERTIONS)[colnames(ORIGRECORDS_INSERTIONS) == "LEN"] <- "EVENT_LENGTH"
      ORIGRECORDS_INSERTIONS$SEPARATOR <- ""
      ORIGRECORDS_INSERTIONS <- ORIGRECORDS_INSERTIONS[,c("CHROMOSOME","START" ,"STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                                          "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                                          "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents",
                                                          "SEPARATOR",
                                                          "EVENT_CALLER", "START_EVENT","END_EVENT","SVTYPE",              
                                                          "VCF_INFO","EVENT_LENGTH" )]
      
      ORIGRECORDS_INSERTIONS <- unique(ORIGRECORDS_INSERTIONS)
      
      #######################################################################################################################
      ##INVERSIONS
      y <- which(INVERSIONS>0)
      START <- y[!(y-1) %in% y]
      STOP <- y[!(y+1) %in% y]
      INVERSIONS_DF <- as.data.frame(cbind(START, STOP))
      INVERSIONS_DF$LEN <- INVERSIONS_DF$STOP - INVERSIONS_DF$START+1
      INVERSIONS_DF$NumOfCallers <- 0
      INVERSIONS_DF <- INVERSIONS_DF[INVERSIONS_DF$LEN>min_sv_length,]
      INVERSIONS_DF$Callers <- ""
      
      INVERSIONS_DF$cnproscanPortion <- 0
      INVERSIONS_DF$delly2Portion <- 0
      INVERSIONS_DF$lumpyPortion <- 0
      INVERSIONS_DF$cnvnatorPortion <- 0
      INVERSIONS_DF$pindelPortion <- 0
      
      INVERSIONS_DF$cnproscanSubevents <- 0
      INVERSIONS_DF$delly2Subevents <- 0
      INVERSIONS_DF$lumpySubevents <- 0
      INVERSIONS_DF$cnvnatorSubevents <- 0
      INVERSIONS_DF$pindelSubevents <- 0
      
      # GET INFO
      for (j in 1:nrow(INVERSIONS_DF)){
        INVERSIONS_DF$NumOfCallers[j] <- max(INVERSIONS[(INVERSIONS_DF$START[j]):(INVERSIONS_DF$STOP[j])])
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
      INVERSIONS_DF <- INVERSIONS_DF[,c("CHROM","START","STOP","LEN","SVTYPE","NumOfCallers","Callers",
                                        "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                        "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents")]
      colnames(INVERSIONS_DF) <- c("CHROMOSOME","START","STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                   "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                   "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents")
      
      ## get back original records foverlaps
      data.table::setDT(INVERSIONS_DF)
      data.table::setkey(INVERSIONS_DF,"CHROMOSOME","START","STOP")
      
      ORIGRECORDS_delly <- data.table::foverlaps(TMPdellyDF, INVERSIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      ORIGRECORDS_pindel <- data.table::foverlaps(TMPpindelDF, INVERSIONS_DF, by.x = c("#CHROM","POS","END"),by.y = c("CHROMOSOME","START","STOP"), mult="all",type="any",nomatch=NULL, which=FALSE)
      
      ORIGRECORDS_INVERSIONS <- rbind(ORIGRECORDS_delly,ORIGRECORDS_pindel,fill=TRUE) #ORIGRECORDS_lumpy
      
      colnames(ORIGRECORDS_INVERSIONS)[colnames(ORIGRECORDS_INVERSIONS) == "#CHROM"] <- "CHROMOSOME"
      colnames(ORIGRECORDS_INVERSIONS)[colnames(ORIGRECORDS_INVERSIONS) == "POS"] <- "START_EVENT"
      colnames(ORIGRECORDS_INVERSIONS)[colnames(ORIGRECORDS_INVERSIONS) == "END"] <- "END_EVENT"
      colnames(ORIGRECORDS_INVERSIONS)[colnames(ORIGRECORDS_INVERSIONS) == "CALLER"] <- "EVENT_CALLER"
      colnames(ORIGRECORDS_INVERSIONS)[colnames(ORIGRECORDS_INVERSIONS) == "LEN"] <- "EVENT_LENGTH"
      ORIGRECORDS_INVERSIONS$SEPARATOR <- ""
      ORIGRECORDS_INVERSIONS <- ORIGRECORDS_INVERSIONS[,c("CHROMOSOME","START" ,"STOP","LENGHT_MERGED","SVTYPE_MERGED","NumOfCallers","Callers",
                                                          "cnproscanPortion","delly2Portion","lumpyPortion","cnvnatorPortion","pindelPortion",
                                                          "cnproscanSubevents","delly2Subevents","lumpySubevents","cnvnatorSubevents","pindelSubevents",
                                                          "SEPARATOR", "EVENT_CALLER", "START_EVENT","END_EVENT","SVTYPE",              
                                                          "VCF_INFO","EVENT_LENGTH" )]
      
      ORIGRECORDS_INVERSIONS <- unique(ORIGRECORDS_INVERSIONS)
      
      # GET TOGETHER # WHAT ABOUT LARGE FILES?
      RESULTS <- rbindlist(list(RESULTS, ORIGRECORDS_DELETIONS,ORIGRECORDS_DUPLICATIONS,ORIGRECORDS_INSERTIONS,ORIGRECORDS_INVERSIONS),fill=TRUE)
      SHORT_RESULTS <- rbindlist(list(SHORT_RESULTS, DELETIONS_DF, DUPLICATIONS_DF, INSERTIONS_DF, INVERSIONS_DF),fill=TRUE)
    }
    
  }
  
  #sort
  RESULTS <- RESULTS[order(RESULTS$CHROMOSOME, RESULTS$START),]
  SHORT_RESULTS <- SHORT_RESULTS[order(SHORT_RESULTS$CHROMOSOME, SHORT_RESULTS$START),]
  
  #write TABLE
  
  write.table(RESULTS, file=output, sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE )
  write.table(SHORT_RESULTS, file=output_short, sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE )
  
}


# develop and test
# args <- c("results/vcf_merged/S01_L001.procaryaSV_merge.tsv",
#           "KP_ref.fasta",
#           "results/lumpy/S01_L001/S01_L001.vcf",
#           "results/delly2/S01_L001/S01_L001.vcf",
#           "results/cnvnator/S01_L001/S01_L001.vcf",
#           "results/pindel/S01_L001/S01_L001.vcf",
#           "results/cnproscan/S01_L001/S01_L001.vcf"
# )
# setwd("/home/rj/1TB/ProcaryoSV_test/")

#run as Rscript

args <- commandArgs(trailingOnly = T)
run_all(args)