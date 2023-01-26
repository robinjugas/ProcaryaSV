#!/bin/bash
# Describe

#conda ENV
eval "$(conda shell.bash hook)"
conda activate aligners

genmap index -F KP_ref.fasta -I mapp_index

genmap map -K 30 -E 2 -I mapp_index -O mapp_genmap -t -w -bg

#  -K, --length INTEGER
#           Length of k-mers
# -E, --errors INTEGER
#           Number of errors


conda deactivate
