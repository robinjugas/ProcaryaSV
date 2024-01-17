
# procaryaSV  MERGE
rule procaryaSV_callers_merge:
    input:
        vcfs=expand("results/{cnv_caller}/{{SAMPLE}}/{{SAMPLE}}.vcf",cnv_caller=SVcallers),
        reference=os.path.join("results/references",config["genome_fasta"])
    output:
        tsv="results/merged_procaryaSV/{SAMPLE}.procaryaSV_callers_merge.tsv",
        venn_png="results/merged_procaryaSV/{SAMPLE}.procaryaSV_venn.png",
        barchart_png="results/merged_procaryaSV/{SAMPLE}.procaryaSV_sv_types.png"
    params:
        sample_name="{SAMPLE}",
        min_sv_length=config["procaryaSV_min_sv_length"],
        max_sv_length=config["procaryaSV_max_sv_length"], # NA means default value which is 1/3 of reference genome length
        minCallersCNV=config["procaryaSV_minCallersCNV"], # for CNV called by 5 tools [DEL,DUP]
        minCallersSV=config["procaryaSV_minCallersSV"], # for SVs called by ~3 tools [INS, INV]
        maxGap=config["procaryaSV_maxGap"],
    log:
        "logs/procaryaSV_callers_merge/{SAMPLE}.log",
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/procaryaSV.yaml")
    script: os.path.join(workflow.basedir,"wrappers/procaryaSV_callers_merge/script.py")

####################################################################################################
# # survivor_merge - remove breseq which generates GD but not VCF
# def survivor_merge_input(wildcards):
#     if "breseq" in SVcallers: #true
#         SVcallersSURVIVOR = SVcallers.remove("breseq")
#         return expand("results/{cnv_caller}/{{SAMPLE}}/{{SAMPLE}}.vcf",cnv_caller=SVcallersSURVIVOR)
#     else: #false
#         return expand("results/{cnv_caller}/{{SAMPLE}}/{{SAMPLE}}.vcf",cnv_caller=SVcallers)

# SURVIVOR MERGE
rule survivor_merge:
    input:
        vcfs=expand("results/{cnv_caller}/{{SAMPLE}}/{{SAMPLE}}.vcf",cnv_caller=SVcallers)
    output:
        vcf="results/merged_survivor/{SAMPLE}.survivor_merged.vcf",
    params:
        sample_files="results/merged_survivor/ls_{SAMPLE}.txt",
        max_allowed_space=config["survivor_max_allowed_space"],#1000,
        min_callers=config["survivor_minCallers"],
        agree_on_type=1,
        agree_on_strand=0,
        estimate_sv_distance=0,
        min_length=config["survivor_min_sv_length"],
    log:
        "logs/survivor_merge/{SAMPLE}.log",
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/survivor.yaml")
    script:
        os.path.join(workflow.basedir,"wrappers/survivor/script.py")


####################################################################################################
# # procaryaSV SAMPLES MERGE
# rule procaryaSV_samples_merge:
#     input:
#         merged_tsvs=expand("results/merged_procaryaSV/{SAMPLE}.procaryaSV_callers_merge.tsv",SAMPLE=config["samples"]),
#         reference=os.path.join("results/references",config["genome_fasta"])
#     output:
#         tsv="results/merged_procaryaSV/{SAMPLE}.procaryaSV_samples_merge.tsv",
#         heatmap_png="results/merged_procaryaSV/{SAMPLE}.procaryaSV_heatmap.png",
#     params:
#         min_overlap=10,
#         max_sv_length="NA",
#         mergeCloseEventsDistanceThreshold=100,
#     log:
#         "logs/procaryaSV_samples_merge/{SAMPLE}.log",
#     threads: 6
#     conda:
#         os.path.join(workflow.basedir,"envs/procaryaSV.yaml")
#     script: os.path.join(workflow.basedir,"wrappers/procaryaSV_samples_merge/script.py")

