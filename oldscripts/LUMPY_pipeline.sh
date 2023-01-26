#!/bin/bash  
# Pipeline for BWA mapping and then LUMPY

#conda
# eval "$(conda shell.bash hook)"
# conda activate quality

REF="KP_ref.fasta"
# #trimming
# for R1 in *R1*.fastq
#     do
#         echo "${R1} trim galore"
#         R2=`echo $R1 | sed 's/_R1_/_R2_/'`
#         bname=`echo $R1 | sed 's/_R1_.\+//'`
#         trim_galore --illumina --quality 20 --paired $R1 $R2
#     done

# #rename
# for filename in *.fq
#     do
#         mv "./$filename" "./$(echo $filename | sed -e 's/_001_val_.//g')";  
#     done 

#conda
# conda deactivate quality

eval "$(conda shell.bash hook)"
conda activate CNVenv

#BWA samtools
bwa index -a is $REF
samtools faidx $REF

mkdir BAMS_folder
mkdir LUMPY_folder

for R1 in *R1.fq
    do
        R2=`echo $R1 | awk -F '_' '{print $1"_L001_R2.fq"}' `
        bname=`echo $R1 | awk -F '_' '{print $1}'`
        #echo $R1
        #echo $R2
        echo "${bname} BWA"

        #LUMPY PREPARATION
        # Align the data
        bwa mem -R "@RG\tID:id\tSM:sample\tLB:lib" -v 1 $REF $R1 $R2 \
            | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
            | samtools view -S -b - \
            > $bname.bam       

        # Extract the discordant paired-end alignments.
        samtools view -b -F 1294 $bname.bam > $bname.discordants.unsorted.bam


        # Extract the split-read alignments
        samtools view -h $bname.bam \
            | lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
            | samtools view -Sb - \
            > $bname.splitters.unsorted.bam

        # Sort both alignments
        samtools sort $bname.discordants.unsorted.bam -o $bname.discordants.bam
        samtools sort $bname.splitters.unsorted.bam -o $bname.splitters.bam
  
        #LUMPY
        lumpyexpress \
            -B $bname.bam \
            -S $bname.splitters.bam \
            -D $bname.discordants.bam \
            -o LUMPY_folder/$bname.vcf

        svtools vcftobedpe -i LUMPY_folder/$bname.vcf -o LUMPY_folder/$bname.bedpe
        #svtools bedpesort LUMPY_folder/$bname.bedpe

        # MOVE TO FINALIZED
        mv *.bam BAMS_folder

        echo "${bname} BWA and LUMPY FINISHED"

    done

#conda
conda deactivate CNVenv





