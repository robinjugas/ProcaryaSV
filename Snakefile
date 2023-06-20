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
include: "rules/breseq.smk"
####################################
# RULE ALL
rule all:
    input:
        #vcf files
        vcf_lumpy=expand("results/lumpy/{SAMPLE}/{SAMPLE}.vcf", SAMPLE = config["samples"]),
        vcf_delly=expand("results/delly2/{SAMPLE}/{SAMPLE}.vcf", SAMPLE = config["samples"]),
        vcf_cnvnator=expand("results/cnvnator/{SAMPLE}/{SAMPLE}.vcf",SAMPLE = config["samples"]),
        vcf_cnproscan=expand("results/cnproscan/{SAMPLE}/{SAMPLE}.vcf",SAMPLE = config["samples"]),
        vcf_pindel=expand("results/pindel/{SAMPLE}/{SAMPLE}.vcf",SAMPLE = config["samples"]),
        # vcf_breseq=expand("results/breseq/{SAMPLE}/output/output.vcf",SAMPLE = config["samples"]),
        # merged SVs
        merge_tsv=expand("results/merged_procaryaSV/{SAMPLE}.procaryaSV_callers_merge.tsv",SAMPLE = config["samples"]),
        vcf_surv=expand("results/merged_survivor/{SAMPLE}.survivor_merged.vcf",SAMPLE = config["samples"]),
        #multiqc
        multiqc="results/reports/ALL_SAMPLES_MULTIQC/ALL_SAMPLES.html"

