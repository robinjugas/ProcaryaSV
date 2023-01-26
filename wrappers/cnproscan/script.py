######################################
# wrapper for rule: CNproScan
######################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## CNproScan \n##\n")
f.close()

shell.executable("/bin/bash")



command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/runCNproScan.R "+\
            snakemake.input.bam + " " +\
            snakemake.input.ref + " " +\
            snakemake.input.coverage + " " +\
            snakemake.input.bedgraph + " " +\
            str(snakemake.threads) + " " +\
            snakemake.output.vcf + " " +\
            snakemake.output.tsv + " " +\
            str(snakemake.params.GCnormalization) + " " +\
            str(snakemake.params.MAPPABILITYnormalization) + " " +\
            str(snakemake.params.ORICnormalization) + " " +\
            str(snakemake.params.ORIC_position) + " " +\
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)
