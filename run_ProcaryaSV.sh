#!/bin/bash
# Run snake file wih set config file, directory and copy files neccessary

WORK_DIR=example/DATA/

WORKFLOW_DIR=example/ProcaryaSV/

CPU_CORES=12

## activate conda ENV
eval "$(conda shell.bash hook)"
conda activate example_snakemake

## RUN SNAKEMAKE

#dry run
snakemake --cores $CPU_CORES -p -n --rulegraph --snakefile $WORKFLOW_DIR/Snakefile --directory $WORK_DIR --configfile config.yaml --dag | dot -Tsvg > dag.svg

# run
snakemake --cores $CPU_CORES --snakefile $WORKFLOW_DIR/Snakefile --directory $WORK_DIR --configfile config.yaml --use-conda --conda-frontend mamba -r

conda deactivate
