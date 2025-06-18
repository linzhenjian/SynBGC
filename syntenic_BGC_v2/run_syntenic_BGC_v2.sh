#!/bin/bash

#SBATCH --account=schmidt-np    # account - abbreviated by -A
#SBATCH --partition=schmidt-np  # partition, abbreviated by -p
export PATH="/uufs/chpc.utah.edu/common/home/schmidt-group3/software/syntenic_BGC_v2":$PATH
#1.make_MCscanx_input.sh -a list -o testout
#1.make_MCscanx_input.sh -i Axi-Age.list -g ./folder_for_ncbi_gbk_files -o test_out
#2.run_interproscan.sh -i test_out
#3.pair_blastp.sh -i test_out
#4.run_mcscanx.sh -i testout 
5.run_bigscape.sh -i testout/ -t strain_taxo
