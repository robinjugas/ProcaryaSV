######################################
# wrapper for rule: procaryaSV_callers_merge
######################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## procaryaSV_callers_merge \n##\n")
f.close()

shell.executable("/bin/bash")


command = "Rscript " + os.path.abspath(os.path.dirname(__file__))+"/procaryaSV_callers_merge.R " +\
    snakemake.output.tsv + " " + \
    snakemake.output.venn_png + " " + \
    snakemake.output.barchart_png + " " + \
    snakemake.input.reference + " " + \
    snakemake.params.sample_name + " " + \
    str(snakemake.params.min_sv_length) + " " + \
    str(snakemake.params.max_sv_length) + " " + \
    str(snakemake.params.minCallersCNV) + " " + \
    str(snakemake.params.minCallersSV) + " " + \
    str(snakemake.params.maxGap) + " " + \
    " ".join(snakemake.input.vcfs) +\
    " >> " + log_filename + " 2>&1"


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()


shell(command)
