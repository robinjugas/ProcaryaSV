#!/bin/bash  
# Skript pro mapování genomu k referenci 

eval "$(conda shell.bash hook)"
conda activate CNVenv


# Germline SV calling
# SV calling is done by sample for high-coverage genomes or in small batches for low-coverage genomes

delly call -g ../reference/FN433596.fasta -o SV.bcf CNVseq.sorted.bam

bcftools view SV.bcf > SV.vcf

conda deactivate


