#!/bin/bash
# Run snake file wih set config file, directory and copy files neccessary

# ## WORK DIR
# WORK_DIR=`pwd`

WORK_DIR=/home/rj/1TB/ProcaryaSV_K.pneumoniae

WORKFLOW_DIR=/home/rj/ownCloud/WORK/ProcaryaSV


## activate conda ENV
eval "$(conda shell.bash hook)"
conda activate snakemake

## RUN SNAKEMAKE

#dry run
snakemake --cores 12 -p -n --rulegraph --snakefile $WORKFLOW_DIR/Snakefile --directory $WORK_DIR --configfile config.yaml --dag | dot -Tsvg > dag.svg

# run
snakemake --cores 12 --snakefile $WORKFLOW_DIR/Snakefile --directory $WORK_DIR --configfile config.yaml --stats $WORK_DIR/snakemake_run.stats.json --use-conda --conda-frontend mamba -r

conda deactivate
