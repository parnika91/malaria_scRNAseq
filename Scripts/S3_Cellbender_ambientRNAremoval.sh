#!/bin/bash/

# conda create -n cellbender python=3.7
# conda activate cellbender
# pip install cellbender

lib=EXP1_l2
cd $lib/h5/

cellbender remove-background \
                 --input raw_feature_bc_matrix.h5 \
                 --output raw_feature_$lib\_out.h5
                 
# From python with PyTables:
# ptrepack --complevel 5 output raw_feature_$lib\_out.h5:/matrix output raw_feature_$lib\_out_seurat.h5:/matrix
