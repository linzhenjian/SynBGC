#!/bin/bash

#SBATCH --account=schmidt-np
#SBATCH --partition=schmidt-np
export PATH="/uufs/chpc.utah.edu/common/home/schmidt-group3/software/MCScanX":$PATH
# --- HELP MESSAGE ---
if [[ -z "$1" || "$1" == "-h" || "$1" == "--help" ]]; then
    cat <<EOF
synteny_BGC, Schmidt Lab, University of Utah
Run MCScanX on paired species comparisons.

Usage: $0 -i <input_directory> 
Note: <input_directory> 
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

#set global viariable for funtion

BIN_PATH="$(cd "$(dirname "$0")" && pwd)"
#BIN_PATH="/scratch/general/vast/zlin/sponge/REF_genome/syntenic_BGC/"
export BIN_PATH
echo "BIN_PATH" $BIN_PATH
echo "input_directory1" $input_directory
input_directory="$(cd "$input_directory" && pwd)"
echo "input_directory2" $input_directory
export input_directory

# --- Check Dependencies ---
for cmd in awk seqkit sed mkdir; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: '$cmd' is required but not found. Aborting." >&2
        exit 1
    fi
done

# --- Prepare Domain File ---
cat "$input_directory"/orign_prot/*_*.tsv > "$input_directory/all.tsv"
awk '{print $1}' "$BIN_PATH/non-BGC-domain" | sed '/^$/d' > "$input_directory/domain_temp"

awk -F '\t' 'NR==FNR{a[$0]; next} !($5 in a)' OFS='\t' "$input_directory/domain_temp" "$input_directory/all.tsv" > "$input_directory/selected.tsv"
awk '{print $1}' "$input_directory/selected.tsv" | awk '!seen[$1]++' > "$input_directory/selected_prot.id"


# --- Pairwise Comparison Loop ---
if [ ! -f "$input_directory/pairs" ]; then
    echo "ERROR: File '$input_directory/pairs' not found." >&2
    exit 1
fi

while read -r name1 name2; do
    [ -z "$name1" ] || [ -z "$name2" ] && continue

    echo "Processing pair: $name1 and $name2"

    gff1="$input_directory/orign_gff/${name1}.gff"
    gff2="$input_directory/orign_gff/${name2}.gff"
    blastp="$input_directory/orign_blast/${name1}_to_${name2}.blastp"

    if [[ ! -f "$gff1" || ! -f "$gff2" || ! -f "$blastp" ]]; then
        echo "WARNING: Missing required files for $name1 and $name2. Skipping."
        continue
    fi

    temp1="$input_directory/${name1}.${name2}.temp1"
    temp2="$input_directory/${name1}.${name2}.temp2"

    cat "$gff1" "$gff2" > "$temp1"
    awk -F '\t' 'NR==FNR{a[$0]}NR>FNR{if ($2 in a) print $0}' "$input_directory/selected_prot.id" "$temp1" > "$input_directory/${name1}_to_${name2}.gff"

    awk -F '\t' 'NR==FNR{a[$0]}NR>FNR{if ($1 in a) print $0}' "$input_directory/selected_prot.id" "$blastp" > "$temp2"
    awk -F '\t' 'NR==FNR{a[$0]}NR>FNR{if ($2 in a) print $0}' "$input_directory/selected_prot.id" "$temp2" > "$input_directory/${name1}_to_${name2}.blast"

    rm -f "$temp1" "$temp2"

    MCScanX "$input_directory/${name1}_to_${name2}" -e 1 -a -s 4 -m 10

done < "$input_directory/pairs"

# --- Combine and Format Output ---
echo "Combining MCScanX output..."

cat "$input_directory"/*.collinearity > "$input_directory/combined"
rm "$input_directory"/*.collinearity "$input_directory"/*.blast "$input_directory"/*.gff

rm -r "$input_directory/collinear"
mkdir -p "$input_directory/collinear"


line_rm=$(cat <<EOF
############### Parameters ###############
# MATCH_SCORE
# MATCH_SIZE
# GAP_PENALTY
# OVERLAP_WINDOW
# E_VALUE
# MAX GAPS
############### Statistics ###############
# Number of collinear genes
# Number of all genes
##########################################
EOF
)

echo "$line_rm" | awk -F ':' 'NR==FNR{a[$0]; next} !($1 in a)' OFS='\t' - "$input_directory/combined" \
    | awk 'BEGIN{n=1} /^## Alignment /{sub(/## Alignment [0-9]+:/,"## Alignment " n ":"); n++} {print}' \
    > "$input_directory/combined.collinear"

# Split alignments into individual files
awk -v outdir="$input_directory/collinear" -F '\t' '
    /^## Alignment/ {part++; next}
    NF >= 3 {
        print $3 > outdir "/part" part ".txt"
        print $2 > outdir "/" part ".txt"
    }
' "$input_directory/combined.collinear"

echo "All done."


#-----------split the syntenic block by gaps >2-----------------

module load parallel  # Load GNU parallel if needed

process_file() {
    file="$1"
    input_dirname=$(dirname "$file")

    # Preprocess input file
    awk -F _ '{print $4, $0}' OFS='\t' "$file" \
        | sed 's/^g//' \
        | sort -k1,1n \
        | awk '{print $2}' \
        | awk '!seen[$0]++' > "$file.tempout"

    mv "$file.tempout" "$file"

    # Extract the name prefix
    name=$(head -n 1 "$file" | awk -F '_' '{print $1"_"$2"_"$3"_"}')
    echo "Processing: $file with name prefix $name"

    # Create the .order file
    grep "$name" "$input_directory/selected_prot.id" \
        | awk -F _ '{print $4, $0}' OFS='\t' \
        | sed 's/^g//' \
        | sort -k1,1n \
        | awk '{print $2}' \
        | awk '!seen[$0]++' > "$file.order"

    # Run the Python script
   python  "$BIN_PATH/script/split_gene_list.py" "$file" "$file.order" "$input_dirname"
}

export -f process_file

# Get number of available CPU cores
cores=$(getconf _NPROCESSORS_ONLN)
echo "Using $cores cores"

# Run *.txt files in parallel
find "$input_directory/collinear" -maxdepth 1 -type f -name '*.txt' -print0 | parallel -0 -j "$cores" process_file


# annotate the genome gbk files
anno_gbk () {
    gbk="$1"
    name=$(basename "$gbk" .gbk)  # correct usage of basename with extension stripping
    awk -F '\t' '{print $1, $5, $6}' OFS='\t' "$input_directory/orign_prot/$name.tsv" \
        | sed 's/_/\t/3' \
        | sed 's/^[^\t]*\t//' \
        | awk '!seen[$0]++' > "$input_directory/$name.temp"
    python "$BIN_PATH/script/annotate_gbk_domain.py" -i "$gbk" -d "$input_directory/$name.temp" -o "$input_directory/${name}_anno.GBK"
    rm  "$input_directory/$name.temp" 
}

# Make sure input_directory, BINPATH, and cores are defined before this
find "$input_directory/" -maxdepth 1 -type f -name '*.gbk' -print0 | parallel -0 -j "$cores" anno_gbk


#-----------add the genes around the syntenic genes-----------------

rm -r "$input_directory/collinear/processed"
mkdir -p "$input_directory/collinear/processed"

add_genes() {
    file="$1"

    # Preprocess input file
    awk -F _ '{print $4,$0}' OFS='\t' "$file" | \
        sed 's/^g//' | \
        sort -k1,1n | \
        awk '{print $2}' | \
        awk '!seen[$0]++' > "$file.tempout"
    mv "$file.tempout" "$file"

    # Extract the name prefix
    name_prefix=$(head -n 1 "$file" | awk -F '_' '{print $1"_"$2"_"$3"_"}')
    echo "Processing: $file with name prefix $name_prefix"

    # Create the .order file
    grep "$name_prefix" "$input_directory/selected_prot.id" | \
        awk -F _ '{print $4,$0}' OFS='\t' | \
        sed 's/^g//' | \
        sort -k1,1n | \
        awk '{print $2}' | \
        awk '!seen[$0]++' > "$file.order"

    # Build safe output path inside processed folder
    input_dirname=$(dirname "$file")
    filename=$(basename "$file")
    output_file="$input_dirname/processed/$filename"

    # Run check script
    OUTPUT=$(python "$BIN_PATH/script/check_neighbour.py" "$file")

    if [[ "$OUTPUT" == *"Valid file"* ]]; then
        echo "Processing file: $file"
        echo "Output file: $output_file"
        python "$BIN_PATH/script/pick_synteny_region.py" "$file" "$file.order" "$output_file"
    fi
}

export -f add_genes

# Run in parallel
find "$input_directory/collinear/" -maxdepth 1 -type f -name '*_part*.txt' -print0 | \
    parallel -0 -j "$cores" add_genes

#------------remove the duplicated blocks in the processed directory---------------

# Find all .txt files in the specified folder (recursively), compute md5sums
find "$input_directory/collinear/processed" -type f -name "*.txt" -exec md5sum {} + > $input_directory/collinear/processed/md5sums

# Process md5sums to detect duplicate checksums
awk '
{
    checksum = $1
    $1 = ""; sub(/^ /, "", $0)
    filename = $0
    if (checksum in arr) {
        arr[checksum] = arr[checksum] ORS filename
    } else {
        arr[checksum] = filename
    }
}
END {
    for (key in arr) {
        num = split(arr[key], files, "\n")
        if (num > 1) {
            print "Duplicate files with checksum " key ":"
            for (i=1; i<=num; i++) print files[i]
            print "------"
        }
    }
}' $input_directory/collinear/processed/md5sums > $input_directory/collinear/processed/duplicated.group

grep -A1 Duplicate $input_directory/collinear/processed/duplicated.group | sed '/Duplicate/d; /--/d; s/ //g'  > $input_directory/collinear/processed/uniquetxt

find  "$input_directory/collinear/processed/" -maxdepth 1  -name "*.txt"  > $input_directory/collinear/processed/all.list
sed -i 's/[.]\///' $input_directory/collinear/processed/all.list

sed  '/^Duplicate/d; /^--/d; s/ //g' $input_directory/collinear/processed/duplicated.group >  $input_directory/collinear/processed/duplicated.list

awk 'NR==FNR{a[$0]=1}NR>FNR{if(a[$0]!=1)print}' $input_directory/collinear/processed/duplicated.list $input_directory/collinear/processed/all.list >> $input_directory/collinear/processed/uniquetxt
awk 'NR==FNR{a[$0]=1}NR>FNR{if(a[$0]!=1)print}' $input_directory/collinear/processed/uniquetxt $input_directory/collinear/processed/all.list > $input_directory/collinear/processed/to_rm

while read i; do rm $i; done < $input_directory/collinear/processed/to_rm

#--------------------add annotation domains to the processed block genes-------------------------


annotate_file() {
    file="$1"
    awk -F '\t' 'NR==FNR{a[$0]}NR>FNR{if ($1 in a) print $0}'  OFS='\t' $file  "$input_directory/selected.tsv"  > "${file}.anno"
}

export -f annotate_file

# Run in parallel
find "$input_directory/collinear/processed" -maxdepth 1 -type f -name '*.txt' -print0 | parallel -0 -j $(nproc) annotate_file


#-------------------------make gbk files from the processed block genes--------------------------------

rm -r  $input_directory/GBK_file
mkdir -p $input_directory/GBK_file
generate_gbk() {
    cluster="$1"
    base=$(basename "$cluster" .txt) #extract the base name of a file, removing the .txt extension

    if [ ! -f "$cluster" ]; then
        echo "Missing file: $cluster"
        return
    fi

    line=$(head -n 1 "$cluster")
    sp_name=$(awk -F _ '{print $1"_"$2}' <<< "$line")
    contig=$(awk -F _ '{print $3}' <<< "$line")
    id1=$(awk -F _ '{print $1"_"$2"_"$3}' <<< "$line")
    id="${id1}_${base}"

    echo "Processing ID: $id"
    python $BIN_PATH/script/extract_contig_gbk.py "$input_directory/orign_gbk/$sp_name.gbk" "$input_directory/$base.gbk" "$contig"
    awk -F '_' '{print $4}' $cluster >  "$input_directory/$base.id"
    python $BIN_PATH/script/extract_records_gbk.py "$input_directory/$base.gbk" "$input_directory/GBK_file/$id.gbk" "$input_directory/$base.id"
    #rm "$input_directory/$base.id" "$input_directory/$base.gbk" 
}

export -f generate_gbk 

find "$input_directory/collinear/processed/" -maxdepth 1 -type f -name '*.txt' -print0 | parallel -0 -j "$cores" generate_gbk






