
import os
import subprocess
import sys
from snakemake.shell import shell

log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: survivor merge \n##\n")
f.close()

shell.executable("/bin/bash")

a_list=snakemake.input.vcfs
textfile = open(snakemake.params.sample_files, "w")
for element in a_list:
    textfile.write(element + "\n")
textfile.close()


command = "SURVIVOR merge " + str(snakemake.params.sample_files) + \
        " " + str(snakemake.params.max_allowed_space) +\
        " " + str(snakemake.params.min_callers) +\
        " " + str(snakemake.params.agree_on_type) +\
        " " + str(snakemake.params.agree_on_strand) +\
        " " + str(snakemake.params.estimate_sv_distance) +\
        " " + str(snakemake.params.min_length) +\
        " " + snakemake.output.vcf +\
        " >> " + log_filename + " 2>&1"


f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)




