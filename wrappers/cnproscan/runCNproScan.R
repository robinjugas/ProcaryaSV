
# SCRIPT TO BE RUN WITH CNPROSCAN v1.0+
run_all <- function(args){
  bamFile <- args[1]
  fastaFile <- args[2]
  coverageFile <- args[3]
  bedgraphFile <- args[4]
  cores <- as.integer(args[5])
  
  VCFfile <- args[6]
  TSVfile <- args[7]
  
  GCnorm <- toupper(args[8])
  MAPnorm <- toupper(args[9])
  ORICnorm <- toupper(args[10])
  oriCposition <- as.integer(args[11])
  
  library(CNproScan)
  DF <- CNproScan::CNproScanCNV(coverageFile, bamFile, fastaFile, GCnorm=GCnorm, MAPnorm=MAPnorm, ORICnorm=ORICnorm,
                     bedgraphFile,oriCposition=, cores=cores)
  # write VCF
  CNproScan::writeVCF(DF,VCFfile)
  # write TSV
  write.table(DF, file = TSVfile, row.names=FALSE, col.names = TRUE, sep="\t")
  
  return(NULL)
}


#run as Rscript
args <- commandArgs(trailingOnly = T)
# install package
if (!require(CNproScan)) {
  devtools::install_github("robinjugas/CNproScan", dependencies = FALSE)
  library(CNproScan)
}


# run main function
run_all(args)
