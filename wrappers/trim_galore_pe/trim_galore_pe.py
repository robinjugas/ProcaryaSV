from snakemake.shell import shell
import os.path
import subprocess
shell.executable("/bin/bash")

log = snakemake.log_fmt_shell()

# Don't run with --fastqc flag
if "--fastqc" in snakemake.params.get("extra", ""):
    raise ValueError("Please remove the --fastqc flag.")

out_dir = os.path.dirname(snakemake.output[0])


command = "trim_galore --fastqc --paired --gzip " + snakemake.params.extra + \
    " " + str(snakemake.input) + " -o " + out_dir + " " + log
shell(command)
