# ProcaryaSV
ProcaryaSV is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to call SV, mainly CNVs from bacterial genomes. 
It employs 5 callers - [DELLY2](https://github.com/dellytools/delly), [LUMPY](https://github.com/arq5x/lumpy-sv),
[Pindel](https://github.com/genome/pindel), [CNVnator](https://github.com/abyzovlab/CNVnator),[INSurVeyor](https://github.com/kensung-lab/INSurVeyor),
and [CNproScan](https://github.com/robinjugas/CNproScan). 

It provides two comparable outputs of merged SV - by ProcaryaSV's merging and by [SURVIVOR's](https://github.com/fritzsedlazeck/SURVIVOR) merge command. 
It starts with optional read trimming, reads alignment and then proceeds to CNV/SV calling. 

## Benchmarking Datasets and Results
All the sequencing reads and scripts are in the Onedrive folder here:

[https://1drv.ms/f/s!Ah7xah3UhCitj4wk3LoPqk9C1OkdxQ?e=IyvG4j](https://1drv.ms/f/s!Ah7xah3UhCitj4wk3LoPqk9C1OkdxQ?e=IyvG4j)

The shared folder contains the "DATASETS_RUN" folder with both dataset, their configuration files, and bash scripts to run them. Just check the bash script and configuration file config.yaml for correct paths leading to FASTA references and ProcaryaSV pipeline folder. 

## Snakemake environment
Snakemake workflow management can run inside [conda](https://docs.conda.io/en/latest/). Create the conda snakemake environemnt easily with already installed conda:
```
conda create --name <SNAKEMAKE_ENVIRONMENT_NAME> -c bioconda -c conda-forge snakemake mamba
```
or setup the environment path:
```
conda create --prefix /<path>/<SNAKEMAKE_ENVIRONMENT_NAME> -c bioconda -c conda-forge snakemake mamba
```

## ProcaryaSV required inputs:
<ul>
<li><strong>yaml configuration file</strong> -- you can find the example in config.yaml with all parameters </li>
<li><strong>sequencing reads</strong> --  place reads in the folder structure: "<EXAMPLE_DATA>/results/raw_reads". The path structure /results/raw_reads is necessary. <EXAMPLE_DATA> is folder name defined by you.  </li>
<li><strong>fasta reference file</strong> -- filename and filepath specified in the yaml configuration file </li> 
</ul>

## Required folder structure
Please, create the initial following structure and place the sequencing data within:
```
<EXAMPLE_DATA>/results/raw_reads    (EXAMPLE_DATA is folder name defined by you)
```

Other folders created later by the pipeline:
```
<EXAMPLE_DATA>/results  (EXAMPLE_DATA is folder name defined by you)
    - raw_reads - created by you, put the sequencing reads here
    - callers_name - created by pipeline - folders for each SV/CNV caller's results
    - references - created by pipeline - alignment indexes
    - mapped - created by pipeline - BAM files
    - merged_procaryaSV - created by pipeline - procaryaSV CNV and SV called
    - merged_survivor - created by pipeline - survivor CNV and SV called
    - reports - created by pipeline - Qualimap, FastQC, Picard reports
<EXAMPLE_DATA>/results/logs - created by pipeline, logs stored here
```

## ProcaryaSV outputs
<ul>
<li> <strong>results/merged_procaryaSV</strong> - ProcaryaSV TSV files with merged SVs feature and coordinates; Venn diagrams; SV types representation plots </li>
<li> <strong>results/merged_survivor</strong> - SURVIVOR VCF files with merged SVs </li>
<li> <strong>results/mapped</strong> - BAM files </li>
<li> <strong>results/reports</strong> - various reports from FASTQC, Picard and Qualimap and </li>
<li> <strong>results/{delly2;lumpy;cnvnator;cnproscan;pindel}</strong> - caller's outputs </li>
</ul>

## How to run the ProcaryaSV pipelie
When you have filled the yaml configuration file and put the files in the folder structure you can start the Snakemake pipeline with this command:
```
snakemake --cores 12 --snakefile path_to_pipeline/Snakefile --directory path_to_DATA --configfile path/config.yaml --use-conda 
```

    --cores 12 - define number of threads
    --snakefile path_to_pipeline/Snakefile - define the path to the ProcaryaSV  Snakefile
    --directory path_to_DATA  - define the path to the data folder structure containing the results/raw_reads
    --configfile path/config.yaml  - define the path to the configuration  Snakefile

    --conda-frontend mamba -r   -  optional arguments, but recommended to use mamba for quicker package installation

Alternatively, you can run the example bash script **run_ProcaryaSV.sh**, which you can modify. 

**If you find some bugs, please fill the github bug report. I'll likely need example of the data to recreate the issue.**

