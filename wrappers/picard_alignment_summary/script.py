######################################
# wrapper for rule: picard_alignment_summary
######################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: picard_alignment_summary by Picard CollectAlignmentSummaryMetrics \n##\n")
f.close()


command = "picard CollectAlignmentSummaryMetrics  " +\
            "I=" + str(snakemake.input.bam) + " " +\
            "R=" + str(snakemake.input.ref) + " " +\
            "O=" + str(snakemake.output.metrics) +\
            " 2>>" + log_filename

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)
