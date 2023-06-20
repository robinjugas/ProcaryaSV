rule samtools_depth:
    input:
        "results/mapped/{SAMPLE}.bam",
    output:
       "results/cnproscan/{SAMPLE}/{SAMPLE}.coverage"
    threads:
        1
    conda:
        os.path.join(workflow.basedir, "envs/samtools.yaml")
    shell:
        """   
        samtools depth -a {input} > {output}
        """

rule genmap:
    input:
        os.path.join("results/references",config["genome_fasta"])
    output:
        mapp_index=directory("results/cnproscan/genmap/mapp_index"),
        mapp_genmap=directory("results/cnproscan/genmap/mapp_genmap"),
        bedgraph="results/cnproscan/genmap/mapp_genmap/" + config["genome_fasta"].split('.fasta')[0] + ".genmap.bedgraph" #or replace config["genome_fasta"].replace(".fasta",".genmap.bedgraph")
    threads:
        1
    log:
        "logs/genmap/genmap.log"
    conda:
        os.path.join(workflow.basedir, "envs/genmap.yaml")
    shell:
        """
        genmap index -F {input} -I {output.mapp_index}  >> {log} 2>&1
        genmap map -K 30 -E 2 -I {output.mapp_index} -O {output.mapp_genmap} -t -w -bg  >> {log} 2>&1
        """


rule cnproscan:
    input:
        bam = "results/mapped/{SAMPLE}.bam",
        ref = os.path.join("results/references",config["genome_fasta"]),
        coverage = "results/cnproscan/{SAMPLE}/{SAMPLE}.coverage",
        bedgraph="results/cnproscan/genmap/mapp_genmap/" + config["genome_fasta"].split('.fasta')[0] + ".genmap.bedgraph"
    output:
        vcf="results/cnproscan/{SAMPLE}/{SAMPLE}.vcf",
        tsv="results/cnproscan/{SAMPLE}/{SAMPLE}.tsv",
    params:
        GCnormalization=config["cnproscan_GC_normalization"],
        MAPPABILITYnormalization=config["cnproscan_MAPPABILITY_normalization"],
        ORICnormalization=config["cnproscan_ORIC_normalization"],
        ORIC_position=config["cnproscan_ORIC_POSITION"]
    threads: 2
    conda:
        os.path.join(workflow.basedir, "envs/cnproscan.yaml")  # add python to ENV
    log:
        "logs/cnproscan/{SAMPLE}_cnproscan.log"
    script: os.path.join(workflow.basedir, "wrappers/cnproscan/script.py")
