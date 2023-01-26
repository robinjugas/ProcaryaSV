##CNOGPRO
library(CNOGpro)
setwd("/home/rj/4TB/SHARED/CNOGpro_artificialdataNEW/umela_data_CNV_coverage200")

# samtools view CNVseq.sorted.bam | perl -lane 'print "$F[2]\t$F[3]"' > out.hits
# samtools view CNVseq.sorted.bam | perl -lane 'print "$F[2]\t$F[3]"' > out.hits

CNOGpro("out.hits", "FN433596.gb", windowlength = 100, name ="FN433596_cov200")


