######################################
# wrapper for rule: picard_gc_bias
######################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: picard_gc_bias by Picard CollectGcBiasMetrics \n##\n")
f.close()


command = "picard CollectGcBiasMetrics  " +\
            "I=" + str(snakemake.input.bam) + " " +\
            "R=" + str(snakemake.input.ref) + " " +\
            "O=" + str(snakemake.output.metrics) + " " +\
            "CHART=" + str(snakemake.output.pdf) + " " +\
            "S=" + str(snakemake.output.summary) +\
            " 2>>" + log_filename

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)
