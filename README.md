# ProcaryaSV
ProcaryaSV is a Snakemake pipeline to call SV, mainly CNVs from bacterial genomes. 
It employs 5 callers - DELLY2, LUMPY, Pindel, CNVnator and CNproScan. 
It provides two comparable outputs of merged SV - by ProcaryaSV's merging and by SURVIVOR's merge command. 


## Inputs:
<ul>
<li><strong>yaml configuration file</strong> -- you can find the example in config.yaml with all parameters </li>
<li><strong>sequencing reads</strong> --  place reads in the folder structure: "EXAMPLE_DATA/results/raw_reads". The path structure /results/raw_reads is necessary. EXAMPLE_DATA is folder name defined by you.  </li>
<li><strong>fasta reference file</strong> -- filename and filepath specified in the yaml configuration file </li> 
</ul>

## EXAMPLE_DATA Folder Structure:
Please, create the initial following structure and place the data within:
```
EXAMPLE_DATA/results/raw_reads    (EXAMPLE_DATA is folder name defined by you)
```

Other folders created lately by pipeline:
```
EXAMPLE_DATA/results  (EXAMPLE_DATA is folder name defined by you)
    - raw_reads - create first and put the sequencing reads here
    - callers_name - created by pipeline - folders for each SV/CNV caller's results
    - references - created by pipeline - alignment indexes
    - mapped - created by pipeline - BAM files
    - merged_procaryaSV - created by pipeline - procaryaSV CNV and SV called
    - merged_survivor - created by pipeline - survivor CNV and SV called
    - reports - created by pipeline - Qualimap, FastQC, Picard reports
--EXAMPLE_DATA/results/logs - created by pipeline, logs stored here
```

## Outputs:
<ul>
<li> <strong>results/merged_procaryaSV</strong> - ProcaryaSV TSV files with merged SVs; Venn diagrams; SV types representation plots </li>
<li> <strong>results/merged_survivor</strong> - SURVIVOR VCF files with merged SVs </li>
<li> <strong>results/mapped</strong> - BAM files </li>
<li> <strong>results/reports</strong> - various reports from FASTQC, Picard and Qualimap and </li>
<li> <strong>results/{delly2;lumpy;cnvnator;cnproscan;pindel}</strong> - caller's outputs </li>
</ul>

## How to RUN the ProcaryaSV pipelie
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

