#!/bin/sh

##SBATCH --account=schmidt-np
##SBATCH --partition=schmidt-np
#SBATCH --ntasks=20
##SBATCH --ntasks-per-node=1  #not use with --mem
#SBATCH --mem=50G

     blastp -query  ./orign_prot/$1.fa -db  ./$2.fa -num_alignments 5  -out  ./orign_blast/${1}_to_${2}.blastp -evalue 1e-10 -num_threads $(nproc --all) -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

