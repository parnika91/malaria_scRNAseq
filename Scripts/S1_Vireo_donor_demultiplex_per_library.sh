#!/bin/bash/

# Vireo run to demultiplex donors beased on SNPs

# install vireo:
# pip install --upgrade --no-deps vireoSNP

# Choose library
lib="EXP1_l1"

# find number of donors in library


# Run vireo with n_donors
# get cellSNP_mat data after running cellSNP on cellranger output
vireo -c cellSNP_mat/ -N $n_donors -o vireo_result/