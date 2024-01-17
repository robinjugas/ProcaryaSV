rule insurveyor:
    input:
        bam = "results/mapped/{SAMPLE}.bam",
        ref = os.path.join("results/references",config["genome_fasta"])
    output:
        vcf="results/insurveyor/{SAMPLE}/{SAMPLE}.vcf",
        dir=directory("results/insurveyor/{SAMPLE}")
    params:
        editedBAM="results/insurveyor/{SAMPLE}/{SAMPLE}.bam"
    threads:
        6
    conda:
        os.path.join(workflow.basedir, "envs/insurveyor.yaml")
    log:
        "logs/insurveyor/{SAMPLE}_insurveyor.log"
    shell:
        """
        picard FixMateInformation I={input.bam} O={params.editedBAM} >> {log} 2>&1
        picard BuildBamIndex I={params.editedBAM} O={params.editedBAM}.bai >> {log} 2>&1
        insurveyor.py --threads {threads} {params.editedBAM} {output.dir} {input.ref} >> {log} 2>&1
        bgzip -d {output.dir}/out.pass.vcf.gz
        cp  {output.dir}/out.pass.vcf {output.vcf}
        """

#https://stackoverflow.com/questions/68378963/snakemake-conda-env-path-and-homer-annotatepeaks-pl
# cp {output.dir}/out.pass.vcf.gz 