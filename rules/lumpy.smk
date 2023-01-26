# Extract the discordant paired-end alignments.
rule samtools_discordant_reads:
    input:
        "results/mapped/{SAMPLE}.bam",
    output:
       "results/lumpy/{SAMPLE}/{SAMPLE}.discordants.bam"
    threads:
        1 
    conda:
        os.path.join(workflow.basedir, "envs/samtools.yaml")
    shell:
        """        
        samtools view -b -F 1294 -@ {threads} -o {wildcards.SAMPLE}.discordants.unsorted.bam {input}
        samtools sort -@ {threads} -T {wildcards.SAMPLE} -o {output} {wildcards.SAMPLE}.discordants.unsorted.bam
        rm {wildcards.SAMPLE}.discordants.unsorted.bam
        """
#
# Extract the split-read alignments
rule samtools_split_reads:
    input:
        "results/mapped/{SAMPLE}.bam",
    output:
        "results/lumpy/{SAMPLE}/{SAMPLE}.splitters.bam"
    params:
        lumpy_script = os.path.join(workflow.basedir,"scripts/lumpy-sv-0.3.1/extractSplitReads_BwaMem")
    threads:
        1
    conda:
        os.path.join(workflow.basedir, "envs/samtools.yaml")
    shell:
        """
        samtools view -h {input} | {params.lumpy_script} -i stdin | samtools view -Sb - > {wildcards.SAMPLE}.splitters.unsorted.bam
        samtools sort -@ {threads} -o {output} {wildcards.SAMPLE}.splitters.unsorted.bam
        rm {wildcards.SAMPLE}.splitters.unsorted.bam
        """

rule lumpy_express:
    input:
        all = "results/mapped/{SAMPLE}.bam",
        discordants = "results/lumpy/{SAMPLE}/{SAMPLE}.discordants.bam",
        splitters = "results/lumpy/{SAMPLE}/{SAMPLE}.splitters.bam"
    output:
        "results/lumpy/{SAMPLE}/{SAMPLE}.vcf"
    threads:
        1
    conda:
        os.path.join(workflow.basedir, "envs/lumpy.yaml")
    log:
        "logs/lumpy/{SAMPLE}_lumpy_express.log"
    shell:
        """
        lumpyexpress -B {input.all} -S {input.splitters} -D {input.discordants} -o {output} >> {log} 2>&1
        """
