import os
import snakemake.io
from multiprocessing import cpu_count


####################################

###########################################
# DEFINITION OF VARIABLES
R1R2 = config["paired_reads_tags"] #["_R1", "_R2"]

SVcallers = config["callers"]

####################################
# SEPARATE RULES
include: "rules/QC_reads_trimming.smk"
include: "rules/multiqc.smk"
include: "rules/alignment.smk"
include: "rules/lumpy.smk"
include: "rules/pindel.smk"
include: "rules/delly.smk"
include: "rules/cnvnator.smk"
include: "rules/cnproscan.smk"
include: "rules/merge_callers.smk"
# include: "rules/breseq.smk"
include: "rules/INSurVeyor.smk"



####################################
# RULE ALL
##### Target rules #####
def get_ruleall_output(wildcards):
    output_list = []
    if config["QCreports"]:
        output_list="results/reports/ALL_SAMPLES_MULTIQC/ALL_SAMPLES.html"
    return output_list

rule all:
    input:
        #vcf files
        expand("results/lumpy/{SAMPLE}/{SAMPLE}.vcf", SAMPLE = config["samples"]),
        expand("results/delly2/{SAMPLE}/{SAMPLE}.vcf", SAMPLE = config["samples"]),
        expand("results/cnvnator/{SAMPLE}/{SAMPLE}.vcf",SAMPLE = config["samples"]),
        expand("results/cnproscan/{SAMPLE}/{SAMPLE}.vcf",SAMPLE = config["samples"]),
        expand("results/pindel/{SAMPLE}/{SAMPLE}.vcf",SAMPLE = config["samples"]),
        expand("results/insurveyor/{SAMPLE}/{SAMPLE}.vcf",SAMPLE = config["samples"]),
        # gd_breseq=expand("results/breseq/{SAMPLE}/{SAMPLE}.gd",SAMPLE = config["samples"]),
        # merged SVs
        expand("results/merged_procaryaSV/{SAMPLE}.procaryaSV_callers_merge.tsv",SAMPLE = config["samples"]),
        expand("results/merged_survivor/{SAMPLE}.survivor_merged.vcf",SAMPLE = config["samples"]),
        #multiqc
        get_ruleall_output

