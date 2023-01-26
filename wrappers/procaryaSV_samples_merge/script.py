######################################
# wrapper for rule: procaryaSV_samples_merge
######################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## procaryaSV_merge \n##\n")
f.close()

shell.executable("/bin/bash")


command = "Rscript " +os.path.abspath(os.path.dirname(__file__))+"/procaryaSV_samples_merge.R "+\
            snakemake.output.tsv + " " + \
            snakemake.input.reference + " " + \
            " ".join(snakemake.input.vcfs) +\
            " >> " + log_filename + " 2>&1"


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()


shell(command)
