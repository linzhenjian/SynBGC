#!/bin/sh

# Display help if no arguments or -h is passed
if [[ -z "$1" || "$1" == "-h" || "$1" == "--help" ]]; then
    cat <<EOF
SynBGC, chmidt Lab, University of Utah
blastp comparison of the protein sequences of the species pairs
Usage: $0  -i <input_directory> #same output_directory as in step 1, making the MCscanx input  files
EOF
    exit
fi


# Parse command line options
while getopts "i:" opt; do
  case $opt in
    i)
      input_directory="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Set default output directory if not provided

# Check for input file
if [ -z "$input_directory" ]; then
    echo "Input directory not specified!"
    exit 1
fi

mkdir -p ${input_directory}/orign_blast

#NOTE: to have best results, compare two stain at a time.
awk '{print $2}' ${input_directory}/pairs | awk '!seen[$1]++' |  while read i; do makeblastdb -in ${input_directory}/orign_prot/$i.fa -dbtype prot -out ${input_directory}/$i.fa; done
while read name1 name2; do
     blastp -query  ${input_directory}/orign_prot/$name1.fa -db  ${input_directory}/$name2.fa -num_alignments 5  -out  ${input_directory}/orign_blast/${name1}_to_${name2}.blastp -evalue 1e-10 -num_threads $(nproc --all) -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
done <  ${input_directory}/pairs

rm ${input_directory}/*_*.pdb ${input_directory}/*_*.phr ${input_directory}/*_*.pin ${input_directory}/*_*.pjs ${input_directory}/*_*.pot ${input_directory}/*_*.psq ${input_directory}/*_*.ptf ${input_directory}/*_*.pto
