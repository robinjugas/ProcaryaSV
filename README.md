# ProcaryaSV
ProcaryaSV is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to call SV, mainly CNVs from bacterial genomes. 
It employs 6 callers - [DELLY2](https://github.com/dellytools/delly), [LUMPY](https://github.com/arq5x/lumpy-sv),
[Pindel](https://github.com/genome/pindel), [CNVnator](https://github.com/abyzovlab/CNVnator),[INSurVeyor](https://github.com/kensung-lab/INSurVeyor),
and [CNproScan](https://github.com/robinjugas/CNproScan). 

It provides two comparable outputs of merged SV - by ProcaryaSV's merging and by [SURVIVOR's](https://github.com/fritzsedlazeck/SURVIVOR) merge command. 
It starts with optional read trimming, reads alignment and then proceeds to CNV/SV calling. 

## Benchmarking Datasets and Results
All the sequencing reads and scripts are in the Onedrive folder here:

[https://1drv.ms/f/s!Ah7xah3UhCitj4wk3LoPqk9C1OkdxQ?e=IyvG4j](https://1drv.ms/f/s!Ah7xah3UhCitj4wk3LoPqk9C1OkdxQ?e=IyvG4j)

The shared folder contains the "DATASETS_RUN" folder with both dataset, their configuration files, and bash scripts to run them. Just check the bash script and configuration file config.yaml for correct paths leading to FASTA references and ProcaryaSV pipeline folder. 

## Wiki
<strong>ProcaryaSV manual can be found here or under Wiki tab above:</strong> \
[https://github.com/robinjugas/ProcaryaSV/wiki/ProcaryaSV-Manual](https://github.com/robinjugas/ProcaryaSV/wiki/ProcaryaSV-Manual)



**If you find some bugs, please fill the github bug report. I'll likely need example of the data to recreate the issue.**

## Citation
ProcaryaSV has been published in BMC Bioinformatics. \
Jugas, R., Vitkova, H. ProcaryaSV: structural variation detection pipeline for bacterial genomes using short-read sequencing. BMC Bioinformatics 25, 233 (2024). https://doi.org/10.1186/s12859-024-05843-1 \
[https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05843-1](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05843-1)
