# ProcaryaSV
ProcaryaSV is a Snakemake pipeline to call SV, mainly CNVs from bacterial genomes. 
It employs 5 callers - DELLY2, LUMPY, Pindel, CNVnator and CNproScan. 
It provides two comparable outputs of merged SV - by ProcaryaSV's merging and by SURVIVOR's merge command. 


## Inputs:
<ul>
<li>yaml configuration file -> example in config.yaml</li>
<li>sequencing reads ->  place them here: EXAMPLE_DATA/results/raw_reads </li>
<li>fasta reference file -> filename and file path specified in config.yaml </li> 
</ul>

## EXAMPLE_DATA Folder Structure:
--EXAMPLE_DATA
    - results - create
        - raw_reads - create path and put the sequencing reads here
    - logs - created by pipeline, logs stored here


## Outputs:
<ul>
<li> results/merged_procaryaSV  -> ProcaryaSV TSV files with merged SVs; Venn diagrams; SV types representation </li>
<li> results/merged_survivor -> SURVIVOR VCF files with merged SVs </li>
<li> results/mapped -> BAM file </li>
<li> results/reports - various reports from FASTQC, Picard and Qualimap and </li>
<li> results/{delly2;lumpy;cnvnator;cnproscan;pindel} -> standalone outputs by SV/CNV detection tools </li>
</ul>

## Running Example
Modify and run the example script file (run_ProcaryaSV.sh) or run the command:
snakemake --cores 12 --snakefile path_to_pipeline/Snakefile --directory path_to_DATA --configfile path/config.yaml --use-conda {--conda-frontend mamba -r} optional arguments, but recommended

