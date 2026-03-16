#!/bin/bash
 

export PATH="/uufs/chpc.utah.edu/common/home/schmidt-group3/software/interproscan/interproscan-5.72-103.0":$PATH
ml jdk/11

# Display help if no arguments or -h is passed
if [[ -z "$1" || "$1" == "-h" || "$1" == "--help" ]]; then
    cat <<EOF
SynBGC, chmidt Lab, University of Utah
annotate the protein sequences by running interproscan
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


find $input_directory/orign_prot -type f -name "*.fa" -print0 | while IFS= read -r -d '' file; do
    echo "Processing: $file"

# Extract the base filename (without the path and file extension).
name=$(basename "$file" | sed 's/\.[^.]*$//')
output_path=$(dirname "$file")

echo "Processing: $name"
echo "Output path: $output_path"

# run InterProScan
interproscan.sh -i "$file" -f tsv -dp -cpu $(nproc --all) -o "$output_path/$name.tsv"

# Post-process the output: sort by gene index.
awk -F '[_\t]' '{print $4, $0}' OFS='\t' "$output_path/$name.tsv" \
  | sed 's/^g//' \
  | sort -k1,1n \
  | sed 's/^[^\t]*\t//' \
  > "$output_path/$name.tsv.temp"

mv "$output_path/$name.tsv.temp" "$output_path/$name.tsv"
done


