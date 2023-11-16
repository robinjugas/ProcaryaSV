# ProcaryaSV
ProcaryaSV is a Snakemake pipeline to call SV, mainly CNVs from bacterial genomes. 
It employs 5 callers - DELLY2, LUMPY, Pindel, CNVnator and CNproScan. 
It provides two comparable outputs of merged SV - by ProcaryaSV's merging and by SURVIVOR's merge command. 


## Inputs:
<ul>
<li>yaml configuration file -> you can find the example in config.yaml with all parameters </li>
<li>sequencing reads ->  place them here: "EXAMPLE_DATA/results/raw_reads" . The path structure /results/raw_reads is necessary.  </li>
<li>fasta reference file -> filename and file path specified in config.yaml </li> 
</ul>

## EXAMPLE_DATA Folder Structure:
Create the following structure and place the data within:
--EXAMPLE_DATA /
    - results /  

        - raw_reads - create first and put the sequencing reads here
        - callers_name - created by pipeline - folders for each SV/CNV caller's results
        - references - created by pipeline - alignment indexes
        - mapped - created by pipeline - BAM files
        - merged_procaryaSV - created by pipeline - procaryaSV CNV and SV called
        - merged_survivor - created by pipeline - survivor CNV and SV called
        - reports - created by pipeline - Qualimap, FastQC, Picard reports

    - logs - created by pipeline, logs stored here


## Outputs:
<ul>
<li> results/merged_procaryaSV  -> ProcaryaSV TSV files with merged SVs; Venn diagrams; SV types representation </li>
<li> results/merged_survivor -> SURVIVOR VCF files with merged SVs </li>
<li> results/mapped -> BAM file </li>
<li> results/reports - various reports from FASTQC, Picard and Qualimap and </li>
<li> results/{delly2;lumpy;cnvnator;cnproscan;pindel} -> standalone outputs by SV/CNV detection tools </li>
</ul>

## How to RUN the ProcaryaSV pipelie
When you have filled the yaml configuration file and put the files in the folder structure you can start the Snakemake pipeline with this command:

snakemake --cores 12 --snakefile path_to_pipeline/Snakefile --directory path_to_DATA --configfile path/config.yaml --use-conda 

--cores 12 - define number of threads
--snakefile path_to_pipeline/Snakefile - define the path to the ProcaryaSV  Snakefile
--directory path_to_DATA  - define the path to the data folder structure containing the results/raw_reads
--configfile path/config.yaml  - define the path to the configuration  Snakefile

{--conda-frontend mamba -r} optional arguments, but recommended to use mamba for quicker package installation

Alternatively, you can run the example bash script run_ProcaryaSV.sh, which you can modify. 

If you find some bugs, please fill the github bug report. I'll likely need example of the data to recreate the issue. 

