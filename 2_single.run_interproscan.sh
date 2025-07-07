#!/bin/bash

##SBATCH --account=schmidt    # account - abbreviated by -A
##SBATCH --partition=kingspeak-shared  # partition, abbreviated by -p
#SBATCH --nodes=1
#SBATCH --ntasks=20
##SBATCH --ntasks-per-node=1  #not use with --mem
#SBATCH --mem=80G
##SBATCH --time=24:00:00
export PATH="/uufs/chpc.utah.edu/common/home/schmidt-group3/software/interproscan/interproscan-5.72-103.0":$PATH
ml jdk/11



prot="$1"

# 提取主文件名（不含路径和扩展名）
name=$(basename "$prot" | sed 's/\.[^.]*$//')
output_path=$(dirname "$prot")

echo "Processing: $name"
echo "Output path: $output_path"

# 运行 InterProScan
interproscan.sh -i "$prot" -f tsv -dp -cpu $(nproc --all) -o "$output_path/$name.tsv"

# 后处理输出：按基因编号排序
awk -F '[_\t]' '{print $4, $0}' OFS='\t' "$output_path/$name.tsv" \
  | sed 's/^g//' \
  | sort -k1,1n \
  | sed 's/^[^\t]*\t//' \
  > "$output_path/$name.tsv.temp"

mv "$output_path/$name.tsv.temp" "$output_path/$name.tsv"



