#!/bin/bash
# Run snake file wih set config file, directory and copy files neccessary

WORK_DIR=example/DATA/

WORKFLOW_DIR=example/ProcaryaSV/

CPU_CORES=12

YAML_CONFIG=config.yaml

## activate conda ENV
eval "$(conda shell.bash hook)"
conda activate ProcaryaSV_conda

## RUN SNAKEMAKE

#dry run
snakemake -p -n --rulegraph --snakefile $WORKFLOW_DIR/Snakefile --directory $WORK_DIR --configfile $YAML_CONFIG --dag | dot -Tsvg > dag.svg

# run
snakemake --cores $CPU_CORES --snakefile $WORKFLOW_DIR/Snakefile --directory $WORK_DIR --configfile $YAML_CONFIG --use-conda --conda-frontend mamba -r

conda deactivate
