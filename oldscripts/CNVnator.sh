#!/bin/bash  
# Skript pro mapování genomu k referenci 

eval "$(conda shell.bash hook)"
conda activate CNVnator

# REF="FN433596.fasta"

# EXTRACTING READ MAPPING FROM BAM/SAM FILES
cnvnator -root sample.root -tree CNVseq.sorted.bam 
# GENERATING A READ DEPTH HISTOGRAM
cnvnator -root sample.root -his 100 -fasta ../reference/FN433596.fasta
# CALCULATING STATISTICS
cnvnator -root sample.root -stat 100
# RD SIGNAL PARTITIONING
cnvnator -root sample.root -partition 100 -ngc
# CNV CALLING
cnvnator -root sample.root -call 100 -ngc > results.cnvnator.txt

# cnvnator2VCF.pl -prefix study1 -reference GRCh37 sample1.cnvnator.out /path/to/individual/fasta_files


conda deactivate