#!/bin/bash  
# Skript pro mapování genomu k referenci 

eval "$(conda shell.bash hook)"
conda activate CNVenv

# samtools faidx ../reference/FN433596.fasta


# "To run Pindel on a BAM-file, use the following command: "
# -i/--config-file
#            the bam config file; either this, a pindel input file, or a pindel 
#            config file is required. Per line: path and file name of bam, insert 
#            size and sample tag.     For example: /data/tumour.bam  400  tumour 

FILE=$(find . -print | grep '\.bam$')
echo "${FILE} 348  test"  > config.txt


pindel -i config.txt -f ../reference/FN433596.fasta -o results -c ALL

pindel2vcf -r ../reference/FN433596.fasta -R FN433596 -d 20100101 -p results_D -v results_Deletions.vcf
pindel2vcf -r ../reference/FN433596.fasta -R FN433596 -d 20100101 -p results_SI -v results_ShortInsertions.vcf
pindel2vcf -r ../reference/FN433596.fasta -R FN433596 -d 20100101 -p results_TD -v results_Tandem_Duplications.vcf
pindel2vcf -r ../reference/FN433596.fasta -R FN433596 -d 20100101 -p results_LI -v results_LongInsertions.vcf


conda deactivate


