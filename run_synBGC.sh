#!/bin/bash

#SBATCH --account=schmidt-np    # account - abbreviated by -A
#SBATCH --partition=schmidt-np  # partition, abbreviated by -p
#SBATCH --mem=0
#SBATCH --ntasks=40
##SBATCH --time=24:00:00
export PATH="/uufs/chpc.utah.edu/common/home/schmidt-group3/software/syntenic_BGC_v2":$PATH
output="./output"
speceis_list="./list"
gbff_path=""
#1.make_MCscanx_input.sh -a "$speceis_list"   -g "$gbff_path"  -o "$output"
#2.prot_annotaion.sh -i "$output"
#3.pair_blastp.sh -i "$output" 
4.run_mcscanx.sh -i "$output" -s 14 -m 50 -x 100
5.run_bigscape.sh -i "$output" 


