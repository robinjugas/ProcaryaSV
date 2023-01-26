#QUALIMAP
rule qualimap:
    input:
        "results/mapped/{SAMPLE}.bam"
    output:
        "results/reports/QUALIMAP/{SAMPLE}_qualimap/qualimapReport.html",
        "results/reports/QUALIMAP/{SAMPLE}_qualimap/genome_results.txt"
    params:
        OUTDIRECTORY="results/reports/QUALIMAP/{SAMPLE}_qualimap"
    threads:
        6
    log:
        "logs/qualimap/{SAMPLE}_qualimap.log"
    conda:
        os.path.join(workflow.basedir, "envs/qualimap.yaml")
    shell:
        "qualimap bamqc -bam {input} -ip -nt {threads} -os -outdir {params.OUTDIRECTORY} --java-mem-size=4G >> {log} 2>&1"

#PICARD from Pindel
rule copy_picard:
    input:
        "results/pindel/{SAMPLE}/{SAMPLE}.insert_size_metrics.txt",
        "results/pindel/{SAMPLE}/{SAMPLE}.alignment_summary.txt"
    output:
        "results/reports/picard/{SAMPLE}.insert_size_metrics.txt",
        "results/reports/picard/{SAMPLE}.alignment_summary.txt"
    threads:
        1
    shell:
        """
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        """

#PICARD GC BIAS report
rule picard_gc_bias:
    input:
        ref=os.path.join("results/references",config["genome_fasta"]),
        bam="results/mapped/{SAMPLE}.bam"
    output:
        metrics="results/reports/picard/{SAMPLE}.gc_bias_metrics.txt",
        pdf="results/reports/picard/{SAMPLE}.gc_bias_metrics.pdf",
        summary="results/reports/picard/{SAMPLE}.summary_metrics.txt"
    log:
        "logs/picard/{SAMPLE}/CollectGcBiasMetrics.log"
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/picard.yaml")
    script:
        os.path.join(workflow.basedir,"wrappers/picard_gc_bias/script.py")


#MULTIQC

# MULTIQC INPUT FUNCTION
def multiqc_input(wildcards):
    if config["reads_trimming"]: #true
        return expand(["results/reports/raw_reads_FASTQC/{SAMPLE}_R1_fastqc.html",
            "results/reports/raw_reads_FASTQC/{SAMPLE}_R2_fastqc.html",
            "results/reports/trimmed_reads_FASTQC/{SAMPLE}_R1_fastqc.html",
            "results/reports/trimmed_reads_FASTQC/{SAMPLE}_R2_fastqc.html",
            "results/reports/QUALIMAP/{SAMPLE}_qualimap/qualimapReport.html",
            "results/reports/picard/{SAMPLE}.insert_size_metrics.txt",
            "results/reports/picard/{SAMPLE}.alignment_summary.txt",
            "results/reports/picard/{SAMPLE}.gc_bias_metrics.txt"], SAMPLE = config["samples"])
    else: #false
        return expand(["results/reports/raw_reads_FASTQC/{SAMPLE}_R1_fastqc.html",
            "results/reports/raw_reads_FASTQC/{SAMPLE}_R2_fastqc.html",            
            "results/reports/QUALIMAP/{SAMPLE}_qualimap/qualimapReport.html",
            "results/reports/picard/{SAMPLE}.insert_size_metrics.txt",
            "results/reports/picard/{SAMPLE}.alignment_summary.txt",
            "results/reports/picard/{SAMPLE}.gc_bias_metrics.txt"], SAMPLE = config["samples"])


rule multiqc_all_reports:
    input: multiqc_input
    output:
        folder=directory("results/reports/ALL_SAMPLES_MULTIQC/"),
        html="results/reports/ALL_SAMPLES_MULTIQC/ALL_SAMPLES.html"
    params:
        basename="ALL_SAMPLES",
        input_folder="results/reports/"
    threads:
        6
    log:
        "logs/multiqc/ALL_SAMPLES_multiqc.log"
    conda:
        os.path.join(workflow.basedir, "envs/multiqc.yaml")
    shell:
        """        
        multiqc --dirs {params.input_folder} --force -o {output.folder} -n {params.basename} >> {log} 2>&1
        """
