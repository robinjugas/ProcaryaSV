#BWA-MEM2
#copy reference file
rule copy_reference:
    input:
        os.path.join(config["reference_folder"],config["genome_fasta"])
    output:
        os.path.join("results/references",config["genome_fasta"])
    shell:
        "cp {input:q} {output:q}"

# Creates a bwa-mem2 index.
rule bwa_mem2_index:
    input:
        os.path.join("results/references",config["genome_fasta"])
    output:
        "results/references/reference_genome.0123",
        "results/references/reference_genome.amb",
        "results/references/reference_genome.ann",
        "results/references/reference_genome.bwt.2bit.64",
        "results/references/reference_genome.pac",
    log:
        "logs/bwa-mem2_index/index_reference.log",
    params:
        prefix="results/references/reference_genome"
    conda:
        os.path.join(workflow.basedir,"envs/bwa.yaml")
    shell:
        """
        bwa-mem2 index -p {params.prefix} {input} >> {log} 2>&1
        """

# BWA INPUT FUNCTION RAW/TRIMMED READS
def bwa_mem_input(wildcards):
    if config["reads_trimming"]: #true
        return ["results/trimmed_reads/{SAMPLE}_R1_trimmed.fastq.gz",
            "results/trimmed_reads/{SAMPLE}_R2_trimmed.fastq.gz"]
    else: #false
        return ["results/raw_reads/{SAMPLE}"+config["paired_reads_tags"][0]+config["reads_extension"],
            "results/raw_reads/{SAMPLE}"+config["paired_reads_tags"][1]+config["reads_extension"]]


# Map reads using bwa-mem2, mark duplicates by samblaster and sort and index by sambamba.
rule bwa_mem:
    input:
        reads=bwa_mem_input,
        idx=expand("results/references/reference_genome.{ext}",ext=["amb", "ann","bwt.2bit.64", "pac"])
    output:
        bam="results/mapped/{SAMPLE}.bam",
        index="results/mapped/{SAMPLE}.bam.bai",
    log:
        "logs/bwa_mem2_sambamba/{SAMPLE}.log",
    params:
        extra="-R '@RG\tID:{SAMPLE}\tSM:{SAMPLE}'", #recommended by LUMPY readme
        sort_extra="-q",  # Extra args for sambamba.
        samblaster_extra="--excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20", #recommended by LUMPY readme
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/bwa.yaml")
    script:
        os.path.join(workflow.basedir, "wrappers/bwa/bwa_mem2.py")
