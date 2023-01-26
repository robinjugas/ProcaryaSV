# FASTQC RAW READS
rule fastqc_before_trimming:
    input:
        "results/raw_reads/{SAMPLE}"+config["paired_reads_tags"][0]+config["reads_extension"],
        "results/raw_reads/{SAMPLE}"+config["paired_reads_tags"][1]+config["reads_extension"]
    output:
        "results/reports/raw_reads_FASTQC/{SAMPLE}"+config["paired_reads_tags"][0]+"_fastqc.html",
        "results/reports/raw_reads_FASTQC/{SAMPLE}"+config["paired_reads_tags"][1]+"_fastqc.html",
    log:
        "logs/fastqc_before_trimming/{SAMPLE}.log"
    threads:
        6
    params:
        FASTQCDIR="results/reports/raw_reads_FASTQC/"
    conda:
        os.path.join(workflow.basedir, "envs/fastqc.yaml")
    shell:
        "fastqc --threads {threads} --quiet {input} -o {params.FASTQCDIR} >> {log} 2>&1"

# RENAME FASTQC BEFORE TRIMMING
rule rename_fastqc_before_trimming:
    input:
        "results/reports/raw_reads_FASTQC/{SAMPLE}"+config["paired_reads_tags"][0]+"_fastqc.html",
        "results/reports/raw_reads_FASTQC/{SAMPLE}"+config["paired_reads_tags"][1]+"_fastqc.html",
    output:
        "results/reports/raw_reads_FASTQC/{SAMPLE}_R1_fastqc.html",
        "results/reports/raw_reads_FASTQC/{SAMPLE}_R2_fastqc.html"
    shell:
        """
        mv {input[0]} {output[0]}
        mv {input[1]} {output[1]}
        """


# TRIMMING RAW READS - ADAPTERS + QUALITY <20
rule trim_galore_pe:
    input:
        ["results/raw_reads/{SAMPLE}"+config["paired_reads_tags"][0]+config["reads_extension"],
        "results/raw_reads/{SAMPLE}"+config["paired_reads_tags"][1]+config["reads_extension"]]
    output:
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][0]+"_val_1.fq.gz",
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][1]+"_val_2.fq.gz",
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][0]+"_val_1_fastqc.html",
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][1]+"_val_2_fastqc.html",
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][0]+"_val_1_fastqc.zip",
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][1]+"_val_2_fastqc.zip"
    threads:
        6
    params:
        extra="--illumina --quality 20 --cores 6"
    log:
        "logs/trim_galore/{SAMPLE}.log"
    conda:
        os.path.join(workflow.basedir, "envs/trim_galore.yaml")
    script:
        os.path.join(workflow.basedir, "wrappers/trim_galore_pe/trim_galore_pe.py")

# RENAME FASTQC OF TRIMMED READS
rule copy_fastqc_reports_after_trimming:
    input:
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][0]+"_val_1_fastqc.html",
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][1]+"_val_2_fastqc.html",
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][0]+"_val_1_fastqc.zip",
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][1]+"_val_2_fastqc.zip"
    output:
        "results/reports/trimmed_reads_FASTQC/{SAMPLE}_R1_fastqc.html",
        "results/reports/trimmed_reads_FASTQC/{SAMPLE}_R2_fastqc.html",
        "results/reports/trimmed_reads_FASTQC/{SAMPLE}_R1_fastqc.zip",
        "results/reports/trimmed_reads_FASTQC/{SAMPLE}_R2_fastqc.zip"
    shell:
        """
        mv {input[0]} {output[0]}
        mv {input[1]} {output[1]}
        mv {input[2]} {output[2]}
        mv {input[3]} {output[3]}
        """

# RENAME TRIMMED READS
rule rename_trimmed_reads:
    input:
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][0]+"_val_1.fq.gz",
        "results/trimmed_reads/{SAMPLE}"+config["paired_reads_tags"][1]+"_val_2.fq.gz"
    output:
        "results/trimmed_reads/{SAMPLE}_R1_trimmed.fastq.gz",
        "results/trimmed_reads/{SAMPLE}_R2_trimmed.fastq.gz"
    shell:
        """
        mv {input[0]} {output[0]}
        mv {input[1]} {output[1]}
        """
