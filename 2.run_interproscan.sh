#!/bin/bash
 
#SBATCH --account=schmidt-np    # account - abbreviated by -A
#SBATCH --partition=schmidt-np  # partition, abbreviated by -p
##SBATCH --ntasks=20



export PATH="/uufs/chpc.utah.edu/common/home/schmidt-group3/software/interproscan/interproscan-5.72-103.0":$PATH
ml jdk/11

# Display help if no arguments or -h is passed
if [[ -z "$1" || "$1" == "-h" || "$1" == "--help" ]]; then
    cat <<EOF
syneny_BGC, chmidt Lab, University of Utah
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

# Check if required tools are available
for cmd in awk seqkit sed mkdir; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "$cmd is required but not found. Aborting."
        exit 1
    }
done


find $input_directory/orign_prot -type f -name "*.fa" -print0 | while IFS= read -r -d '' file; do
    echo "Processing: $file"

# 提取主文件名（不含路径和扩展名）
name=$(basename "$file" | sed 's/\.[^.]*$//')
output_path=$(dirname "$file")

echo "Processing: $name"
echo "Output path: $output_path"

# 运行 InterProScan
interproscan.sh -i "$file" -f tsv -dp -cpu $(nproc --all) -o "$output_path/$name.tsv"

# 后处理输出：按基因编号排序
awk -F '[_\t]' '{print $4, $0}' OFS='\t' "$output_path/$name.tsv" \
  | sed 's/^g//' \
  | sort -k1,1n \
  | sed 's/^[^\t]*\t//' \
  > "$output_path/$name.tsv.temp"

mv "$output_path/$name.tsv.temp" "$output_path/$name.tsv"
done


