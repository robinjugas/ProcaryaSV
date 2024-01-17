# breseq
rule breseq:
    input:
        reference=os.path.join("results/references",config["genome_fasta"]),
        reads=bwa_mem_input
    output:
        vcf="results/breseq/{SAMPLE}/{SAMPLE}.vcf",
        gd="results/breseq/{SAMPLE}/{SAMPLE}.gd",
        folder=directory("results/breseq/{SAMPLE}/")
    log:
        "logs/breseq/{SAMPLE}_breseq.log"
    params:
        prefix="{SAMPLE}",
        extra="--cnv",
        tmpOutputVCF="results/breseq/{SAMPLE}/output/output.vcf",
        tmpOutputGD="results/breseq/{SAMPLE}/output/output.gd"
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/breseq.yaml")
    shell:
        """
        breseq -r {input.reference} {input.reads} -j {threads} -n {params.prefix} -o {output.folder} {params.extra} >> {log} 2>&1
        cp {params.tmpOutputGD} {output.gd}
        cp {params.tmpOutputGD} {output.vcf}
        """
