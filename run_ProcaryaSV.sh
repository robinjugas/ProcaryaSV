#!/bin/bash
# Run snake file wih set config file, directory and copy files neccessary

WORK_DIR=example/DATA/

WORKFLOW_DIR=example/ProcaryaSV/


## activate conda ENV
eval "$(conda shell.bash hook)"
conda activate example_snakemake

## RUN SNAKEMAKE

#dry run
snakemake --cores 12 -p -n --rulegraph --snakefile $WORKFLOW_DIR/Snakefile --directory $WORK_DIR --configfile config.yaml --dag | dot -Tsvg > dag.svg

# run
snakemake --cores 12 --snakefile $WORKFLOW_DIR/Snakefile --directory $WORK_DIR --configfile config.yaml --use-conda --conda-frontend mamba -r

conda deactivate
