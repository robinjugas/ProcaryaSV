# FASTQC RAW READS
rule fastqc_before_trimming:
    input:
        "results/raw_reads/{SAMPLE}_R1.fastq.gz",
        "results/raw_reads/{SAMPLE}_R2.fastq.gz"
    output:
        "results/reports/{SAMPLE}/{SAMPLE}_R1_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_R2_fastqc.html"
    log:
        "logs/fastqc_before_trimming/{SAMPLE}.log"
    threads:
        4
    params:
        FASTQCDIR="results/reports/{SAMPLE}"
    conda:
        os.path.join(workflow.basedir, "envs/fastqc.yaml")
    shell:
        "fastqc --threads {threads} --quiet {input} -o {params.FASTQCDIR} >> {log} 2>&1"

# TRIMMING RAW READS - ADAPTERS + QUALITY <20
rule trim_galore_pe:
    input:
        "results/raw_reads/{SAMPLE}_R1.fastq.gz",
        "results/raw_reads/{SAMPLE}_R2.fastq.gz"
    output:
        "results/trimmed_reads/{SAMPLE}_R1_val_1.fq.gz",
        "results/trimmed_reads/{SAMPLE}_R1.fastq.gz_trimming_report.txt",
        "results/trimmed_reads/{SAMPLE}_R2_val_2.fq.gz",
        "results/trimmed_reads/{SAMPLE}_R2.fastq.gz_trimming_report.txt"
    threads:
        2
    params:
        extra="--illumina --quality 20 --cores 2"
    log:
        "logs/trim_galore/{SAMPLE}.log"
    conda:
        os.path.join(workflow.basedir, "envs/trim_galore.yaml")
    script:
        os.path.join(workflow.basedir, "scripts/trim_galore_pe.py")

# FASTQC TRIMMED READS
rule fastqc_after_trimming:
    input:
        "results/trimmed_reads/{SAMPLE}_R1_val_1.fq.gz",
        "results/trimmed_reads/{SAMPLE}_R2_val_2.fq.gz"
    output:
        "results/reports/{SAMPLE}/{SAMPLE}_R1_val_1_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_R2_val_2_fastqc.html"
    log:
        "logs/fastqc_after_trimming/{SAMPLE}.log"
    threads:
        4
    params:
        FASTQCDIR="results/reports/{SAMPLE}"
    conda:
        os.path.join(workflow.basedir, "envs/fastqc.yaml")
    shell:
        "fastqc --threads {threads} --quiet {input} -o {params.FASTQCDIR} >> {log} 2>&1"



rule multiqc:
    input:
        "results/reports/{SAMPLE}/{SAMPLE}_R1_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_R2_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_R1_val_1_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_R2_val_2_fastqc.html",
        "results/reports/bismark_not_deduplicated/{SAMPLE}/{SAMPLE}.bismark2report.html",
        "results/reports/qualimap/{SAMPLE}/qualimapReport.html"
    output:
        directory("results/reports/multiqc/{SAMPLE}")
    params:
        basename="{SAMPLE}",
        INDIRECTORY="results/reports/bismark_not_deduplicated/{SAMPLE}/"
        # DIRECTORY1="results/reports/bismark_deduplicated/{SAMPLE}",
    threads:
        6
    log:
        "logs/multiqc/{SAMPLE}_multiqc.log"
    conda:
        os.path.join(workflow.basedir, "envs/multiqc.yaml")
    shell:
        """
        mkdir {output}
        multiqc --force -o {output} -n {params.basename} {params.INDIRECTORY} >> {log} 2>&1
        """
