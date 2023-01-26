rule delly2:
    input:
        bam = "results/mapped/{SAMPLE}.bam",
        ref = os.path.join("results/references",config["genome_fasta"])
    output:
        vcf="results/delly2/{SAMPLE}/{SAMPLE}.vcf",
        bcf="results/delly2/{SAMPLE}/{SAMPLE}.bcf",
    threads:
        1
    conda:
        os.path.join(workflow.basedir, "envs/delly.yaml")
    log:
        "logs/delly2/{SAMPLE}_delly.log"
    shell:
        """
        delly call -o {output.bcf} -g {input.ref} {input.bam} >> {log} 2>&1
        bcftools view {output.bcf} > {output.vcf}    
        """