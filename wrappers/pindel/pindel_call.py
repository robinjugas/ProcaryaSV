
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: pindel \n##\n")
f.close()

command = "pindel" +\
          " -T " + str(snakemake.threads) +\
          " -i " + snakemake.input.config +\
          " -f " + snakemake.input.ref +\
          " -o " + snakemake.params.prefix +\
          " &> " + log_filename

f = open(log_filename, 'at')
f.write("## COMMAND: " +command+ "\n")
f.close()

shell(command)
