#!/bin/bash


export PATH="/uufs/chpc.utah.edu/common/home/schmidt-group3/software/BiG-SCAPE-1.1.5":$PATH
source activate /uufs/chpc.utah.edu/common/home/schmidt-group2/software/miniconda3/envs/bigscape
BIN_PATH="$(cd "$(dirname "$0")" && pwd)"
# --- HELP MESSAGE ---
if [[ -z "$1" || "$1" == "-h" || "$1" == "--help" ]]; then
    cat <<EOF
SynBGC, Schmidt Lab, University of Utah
Run bigscape  on gbk files.

Usage: $0 -i <input_directory> 
Note: <input_directory>  should be the same ones where MCScanX input files were created.
EOF
    exit 0
fi

# --- Parse Arguments ---
while getopts "i:" opt; do
  case $opt in
    i) input_directory="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
  esac
done

# --- Validate Input Directory ---
if [ -z "$input_directory" ]; then
    echo "ERROR: Input directory not specified." >&2
    exit 1
fi

if [ ! -d "$input_directory" ]; then
    echo "ERROR: Directory '$input_directory' does not exist." >&2
    exit 1
fi



#------------------run bigscale------------------------
rm -r $input_directory/output_mibig
bigscape.py  -i $input_directory/GBK_file/ -o $input_directory/output_mibig  --include_singletons --include_gbk_str "_"  --mix --no_classify --cutoffs 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85  --mode global 

 
 # ------------------------------process network files--------------------------
rm -r $input_directory/network
mkdir -p $input_directory/network

# =======================
# Step 1: Collect network files
# =======================
echo "[Step 1] Collecting network files..."
find "$input_directory/output_mibig" -name 'mix*.network' > "$input_directory/network/networkfile"

# =======================
# Step 2: Extract specific network file
# =======================
mix_mini=$(grep "mix_c0.10.network" "$input_directory/network/networkfile")
if [[ ! -f "$mix_mini" ]]; then
  echo "Error: $mix_mini not found."
  exit 1
fi

# =======================
# Step 3: Create edge list
# =======================
echo "[Step 2] Creating edge list..."
awk '{print $1, $2}' "$mix_mini" > "$input_directory/network/net"

# =======================
# Step 4: Run connected component detection
# =======================
echo "[Step 3] Running connected_component.py..."
python $BIN_PATH/script/connected_component.py "$input_directory/network/net" "$input_directory/network/connected_component.txt"
sed -i '/Clustername/d' $input_directory/network/connected_component.txt

# =======================
# Step 5: Group BGCs
# =======================
echo "[Step 5] Grouping BGCs..."
python $BIN_PATH/script/group_BGC.py -i "$input_directory/network/connected_component.txt" -o "$input_directory/network/gene_group"

# =======================
# Step 6: Generate mapping
# =======================
echo "[Step 6] Generating mapping..."
awk '{print $2}' "$input_directory/network/gene_group" > "$input_directory/network/mapping"

# =======================
# Step 7: Copy and rename network files
# =======================
echo "[Step 7] Copying and renaming network files..."
while read path; do
  name=$(basename "$path")
  cp "$path" "$input_directory/network/"
  python $BIN_PATH/script/rename.py -m "$input_directory/network/mapping" -i "$path" -o "$input_directory/network/$name.rename"


  #remove duplicated edges

  awk 'NR==1' "$input_directory/network/$name.rename"  > "$input_directory/network/$name.temp"
  awk 'NR > 1' "$input_directory/network/$name.rename" | sort -t $'\t' -k1,1 -k2,2 -k3,3n  | awk -F '\t' '!seen[$1,$2]++'  >>  "$input_directory/network/$name.temp"
  mv "$input_directory/network/$name.temp" "$input_directory/network/$name.rename"
done < "$input_directory/network/networkfile"
# =======================
# Step 8: Prepare NODE_info from annotations
# =======================
echo "[Step 8] Preparing NODE_info from annotations..."
ANNOT_FILES=$(find "$input_directory/output_mibig" -name 'Network_Annotations_Full.tsv')

awk '{print $1}' $ANNOT_FILES \
  | awk -F _ '{print $1"_"$2, $0}' OFS='\t' \
  | sed '1d' \
  | awk -F '\t' 'NR==FNR{a[$5]=$0;next} NR>FNR{if (a[$1]=="") {print $0,"-";} else {print $0, a[$1];}}' OFS='\t' "$input_directory/strain_taxo" - \
  | awk '{$1=""; print $0}' \
  | sed 's/ //' \
  | sed 's/ /\t/g' > "$input_directory/network/NODE_info"

HEADER="NODE\tclass\torder\tfamily\tgenus\tspecies"
{ echo -e "$HEADER"; cat "$input_directory/network/NODE_info"; } > "$input_directory/network/NODE_info.tmp" && mv "$input_directory/network/NODE_info.tmp" "$input_directory/network/NODE_info"

# =======================
# Step 9: Prepare NODE_rename_info
# =======================
echo "[Step 9] Preparing NODE_rename_info..."
awk -F , '{print $1, $0}' OFS='\t' "$input_directory/network/mapping" > "$input_directory/network/add"
awk -F _ '{print $1"_"$2, $0}' OFS='\t' "$input_directory/network/add" > "$input_directory/network/add_formatted"

HEADER="NODE\tgroup\tclass\torder\tfamily\tgenus\tspecies"
awk -F '\t' 'NR==FNR{a[$5]=$0;next} NR>FNR{if (a[$1]=="") {print $0,"-";} else {print $0, a[$1];}}' OFS='\t' "$input_directory/strain_taxo" "$input_directory/network/add_formatted" \
  | sed 's/^[^\t]*\t//' > "$input_directory/network/NODE_rename_info"

{ echo -e "$HEADER"; cat "$input_directory/network/NODE_rename_info"; } > "$input_directory/network/NODE_rename_info.tmp" && mv "$input_directory/network/NODE_rename_info.tmp" "$input_directory/network/NODE_rename_info"

# =======================
# Final cleanup
# =======================
echo "[Done] All steps completed. Output is in the $input_directory/network directory."
#rm $input_directory/network/connected_component.txt




