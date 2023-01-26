# procaryaSV  CALLERS OVERLAP
rule procaryaSV_callers_merge:
    input:
        vcfs=expand("results/{cnv_caller}/{{SAMPLE}}/{{SAMPLE}}.vcf",cnv_caller=SVcallers),
        reference=os.path.join("results/references",config["genome_fasta"])
    output:
        tsv="results/merged_procaryaSV/{SAMPLE}.procaryaSV_callers_merge.tsv",
        venn_png="results/merged_procaryaSV/{SAMPLE}.procaryaSV_venn.png",
        barchart_png="results/merged_procaryaSV/{SAMPLE}.procaryaSV_sv_types.png"
    params:
        min_sv_length=10,
        max_sv_length="NA",
        mergeCloseEventsDistanceThreshold=100,
    log:
        "logs/procaryaSV_callers_merge/{SAMPLE}.log",
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/procaryaSV.yaml")
    script: os.path.join(workflow.basedir,"wrappers/procaryaSV_callers_merge/script.py")

# procaryaSV SAMPLES MERGE
rule procaryaSV_samples_merge:
    input:
        merged_tsvs=expand("results/merged_procaryaSV/{SAMPLE}.procaryaSV_callers_merge.tsv",SAMPLE=config["samples"]),
        reference=os.path.join("results/references",config["genome_fasta"])
    output:
        tsv="results/merged_procaryaSV/{SAMPLE}.procaryaSV_samples_merge.tsv",
        heatmap_png="results/merged_procaryaSV/{SAMPLE}.procaryaSV_heatmap.png",
    params:
        min_overlap=10,
        max_sv_length="NA",
        mergeCloseEventsDistanceThreshold=100,
    log:
        "logs/procaryaSV_samples_merge/{SAMPLE}.log",
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/procaryaSV.yaml")
    script: os.path.join(workflow.basedir,"wrappers/procaryaSV_samples_merge/script.py")


# SURVIVOR MERGE
rule survivor_merge:
    input:
        vcfs=expand("results/{cnv_caller}/{{SAMPLE}}/{{SAMPLE}}.vcf",cnv_caller=SVcallers),
    output:
        vcf="results/merged_survivor/{SAMPLE}.survivor_merged.vcf",
    params:
        sample_files="results/merged_survivor/ls_{SAMPLE}.txt",
        max_allowed_space=1000,
        min_callers=2,
        agree_on_type=1,
        agree_on_strand=0,
        estimate_sv_distance=0,
        min_length=30,
    log:
        "logs/survivor_merge/{SAMPLE}.log",
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/survivor.yaml")
    script:
        os.path.join(workflow.basedir,"wrappers/survivor/script.py")