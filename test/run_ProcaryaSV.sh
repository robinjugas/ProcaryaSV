#!/bin/bash
# Run snake file wih set config file, directory and copy files neccessary

# ## WORK DIR
WORK_DIR=`pwd`

# WORK_DIR=/home/rj/4TB/PHD_2024/B_RUN_PROCARYASV/CNOGpro_dataset_RUN

WORKFLOW_DIR=/home/rj/1TB/GIT/ProcaryaSV


## activate conda ENV
eval "$(conda shell.bash hook)"
conda activate snakemake

## RUN SNAKEMAKE

#dry run
snakemake --rulegraph --snakefile $WORKFLOW_DIR/Snakefile --directory $WORK_DIR --configfile config.yaml  --dag | dot -Tpdf -o dag_rulegraph.pdf

# run
snakemake --cores 12 --snakefile $WORKFLOW_DIR/Snakefile --directory $WORK_DIR --configfile config.yaml --conda-prefix /home/rj/4TB/CEITEC/snakemake_conda --use-conda




conda deactivate
