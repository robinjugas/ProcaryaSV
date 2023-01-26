
rule picard_insert_size:
    input:
        bam = "results/mapped/{SAMPLE}.bam",
    output:
        metrics="results/pindel/{SAMPLE}/{SAMPLE}.insert_size_metrics.txt",
        histogram="results/pindel/{SAMPLE}/{SAMPLE}.insert_size_histogram.pdf",
    log: "logs/picard/{SAMPLE}/CollectInsertSizeMetrics.log"
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/picard.yaml")
    script:
        os.path.join(workflow.basedir,"wrappers/picard_insert_size/script.py")

rule picard_alignment_summary:
    input:
        ref=os.path.join("results/references",config["genome_fasta"]),
        bam = "results/mapped/{SAMPLE}.bam",
    output:
        metrics="results/pindel/{SAMPLE}/{SAMPLE}.alignment_summary.txt",
    log:
        "logs/picard/{SAMPLE}/CollectAlignmentSummaryMetrics.log"
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/picard.yaml")
    script:
        os.path.join(workflow.basedir,"wrappers/picard_alignment_summary/script.py")
        

rule generate_pindel_config:
    input:
        bam="results/mapped/{SAMPLE}.bam",
        metrics_insert_size="results/pindel/{SAMPLE}/{SAMPLE}.insert_size_metrics.txt",
        metrics_alignment="results/pindel/{SAMPLE}/{SAMPLE}.alignment_summary.txt",
    output:
        config="results/pindel/{SAMPLE}/{SAMPLE}.cfg",
    log:
        "logs/pindel/{SAMPLE}_generate_pindel_config.log",
    threads: 1
    conda:
        os.path.join(workflow.basedir,"envs/pindel.yaml")
    script:
        os.path.join(workflow.basedir,"wrappers/pindel/generate_pindel_config.py")

rule pindel:
    input:
        config="results/pindel/{SAMPLE}/{SAMPLE}.cfg",
        bam="results/mapped/{SAMPLE}.bam",
        bai="results/mapped/{SAMPLE}.bam.bai",
        ref=os.path.join("results/references",config["genome_fasta"]),
    output:
        pindel= expand("results/pindel/{{SAMPLE}}/{{SAMPLE}}_{ext}",ext=["BP","CloseEndMapped","D","INT_final","INV","LI","RP","SI","TD",],),
    params:
        prefix=lambda wildcards: "results/pindel/%s/%s" % (wildcards.SAMPLE,wildcards.SAMPLE),
    log:
        "logs/pindel/{SAMPLE}_pindel.log",
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/pindel.yaml")
    script:
        os.path.join(workflow.basedir,"wrappers/pindel/pindel_call.py")

rule pindel2vcf:
    input:
        pindel=expand("results/pindel/{{SAMPLE}}/{{SAMPLE}}_{ext}",ext=["BP","CloseEndMapped","D","INT_final","INV","LI","RP","SI","TD",],),
        ref=os.path.join("results/references",config["genome_fasta"]),
    output:
        vcf="results/pindel/{SAMPLE}/{SAMPLE}.vcf",
    params:
        refname=config["genome_fasta"].split('.')[0],
        refdate=config.get("pindel_vcf_refdate", "202202"),
    log:
        "logs/pindel/{SAMPLE}_pindel2vcf.log",
    threads: 1
    conda:
        os.path.join(workflow.basedir,"envs/pindel.yaml")
    script:
        os.path.join(workflow.basedir,"wrappers/pindel/pindel2vcf.py")