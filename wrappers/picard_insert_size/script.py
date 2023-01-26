######################################
# wrapper for rule: picard_insert_size
######################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: picard_insert_size by Picard CollectInsertSizeMetrics \n##\n")
f.close()

# java -jar /home/rj/4TB/RELAPS_TREE/picard.jar MergeVcfs I=20wes/modified/20wes.vcf.gz I=36wes/modified/36wes.vcf.gz I=68wes/modified/68wes.vcf.gz O=AB1561.MERGED.vcf
# https://github.com/bioconda/bioconda-recipes/issues/14986

command = "picard CollectInsertSizeMetrics  " +\
            "I=" + str(snakemake.input.bam) + " " +\
            "O=" + str(snakemake.output.metrics) + " " +\
            "H=" + str(snakemake.output.histogram) + " " +\
            "M=0.5"  +\
            " 2>>" + log_filename

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)
