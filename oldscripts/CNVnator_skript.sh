#!/bin/bash
# Skript pro mapování genomu k referenci

eval "$(conda shell.bash hook)"
conda activate CNVnator

REF="KP_ref.fasta"

mkdir VCF_results
mkdir CNVnator_results

for BAM in *.bam; do

    SAMPLENAME=$(echo $BAM | awk -F '.' '{print $1}')

    # EXTRACTING READ MAPPING FROM BAM/SAM FILES
    cnvnator -root sample.root -tree $BAM
    # GENERATING A READ DEPTH HISTOGRAM
    cnvnator -root sample.root -his 100 -fasta $REF
    # CALCULATING STATISTICS
    cnvnator -root sample.root -stat 100
    # RD SIGNAL PARTITIONING
    cnvnator -root sample.root -partition 100 -ngc
    # CNV CALLING
    cnvnator -root sample.root -call 100 -ngc > "CNVnator_results/${SAMPLENAME}_results.txt"

    CNVnator-master/cnvnator2VCF.pl -prefix $SAMPLENAME -reference $REF "CNVnator_results/${SAMPLENAME}_results.txt" /home/rj/4TB/SHARED/CNVNATOR_REAL > "VCF_results/${SAMPLENAME}_results.vcf"

    rm sample.root
    echo "${SAMPLENAME} FINISHED"

done

conda deactivate
