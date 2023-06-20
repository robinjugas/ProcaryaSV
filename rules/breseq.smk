# breseq
rule breseq:
    input:
        reference=os.path.join("results/references",config["genome_fasta"]),
        reads=bwa_mem_input
    output:
        vcf="results/breseq/{SAMPLE}/output/output.vcf",
        folder=directory("results/breseq/{SAMPLE}/")
    log:
        "logs/breseq/{SAMPLE}_breseq.log"
    params:
        prefix="{SAMPLE}",
        extra="--cnv"
    threads: 6
    conda:
        os.path.join(workflow.basedir,"envs/breseq.yaml")
    shell:
        """
        breseq -r {input.reference} {input.reads} -j {threads} -n {params.prefix} -o {output.folder} {params.extra} >> {log} 2>&1
        """