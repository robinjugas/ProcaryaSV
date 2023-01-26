rule cnvnator:
    input:
        bam = "results/mapped/{SAMPLE}.bam",
        ref = os.path.join("results/references",config["genome_fasta"])
    output:
        vcf="results/cnvnator/{SAMPLE}/{SAMPLE}.vcf",
        txt="results/cnvnator/{SAMPLE}/{SAMPLE}.txt",
        root="results/cnvnator/{SAMPLE}/{SAMPLE}.root",
    params:
        bin_size=100,
        reference=config["genome_fasta"],
        ref_folder=config["reference_folder"]
    threads:
        1
    conda:
        os.path.join(workflow.basedir, "envs/cnvnator.yaml")
    log:
        "logs/cnvnator/{SAMPLE}_cnvnator.log"
    shell:
        """        
        cnvnator -root {output.root} -tree {input.bam} >> {log} 2>&1
        cnvnator -root {output.root} -his {params.bin_size} -fasta {input.ref} >> {log} 2>&1
        cnvnator -root {output.root} -stat {params.bin_size} >> {log} 2>&1
        cnvnator -root {output.root} -partition {params.bin_size} -ngc >> {log} 2>&1
        cnvnator -root {output.root} -call {params.bin_size} -ngc > {output.txt} 2>>{log}
        cnvnator2VCF.pl -prefix {wildcards.SAMPLE} -reference {params.reference} {output.txt} {params.ref_folder} > {output.vcf} 2>>{log}
        """

