#!/bin/bash
#SBATCH --account=schmidt-np   # SLURM account name
#SBATCH --partition=schmidt-np # SLURM partition to use

set -euo pipefail              # Safe mode: exit on error, undefined variable, or pipe failure

# Define input and output directories/files
genome_info="genome_info"              # Tab-delimited file: genome path info per strain
GBK="./GBK_for_contig"                 # Output directory for GenBank files
strain_list="Axi-Age.list"             # List of strains and corresponding GFF files
block_dir="collinear/processed/"       # Directory with collinear block files
node_list="$1"                         # Input list of node IDs to process
interpro_anno="./orign_prot"           # Directory with InterProScan annotation files

# Prepare output directory and initialize summary file
> "$GBK/gbkinfo.txt"
mkdir -p "$GBK"

# Main loop over node IDs
while read -r i; do
    # Parse components from the node ID string
    cluster=$(echo "$i" | sed 's/_/ /3' | awk '{print $NF".txt"}')        # Extract cluster name
    mcscan_contig=$(echo "$i" | sed 's/_/ /3' | awk '{print $1}')         # Extract MCSCAN contig
    contig=$(echo "$i" | awk -F '_' '{print $3}')                         # Extract contig ID
    strain_name=$(echo "$i" | awk -F '_' '{print $1"_"$2}')               # Extract strain name

    echo "Processing ID: $i"
    echo "$i" >> "$GBK/gbkinfo.txt"

    # Append collinear block content if file exists
    if [ -f "$block_dir/$cluster" ]; then
        cat "$block_dir/$cluster" >> "$GBK/gbkinfo.txt"
    else
        echo "Missing file: $block_dir/$cluster" >&2
        continue  # Skip to next if missing
    fi

    # Find genome FASTA path and original GFF file path
    genome=$(awk -v name="$strain_name" '$3 == name {print $1}' "$genome_info")
    orign_gff=$(awk -v name="$strain_name" '$1 == name {print $2}' "$strain_list")

    # Create temporary filenames
    temp_prefix="tmp_${i}_$$"                        # Unique prefix using PID
    temp_gff_file="${temp_prefix}_gff.gff"           # Temp GFF for current contig
    temp_anno_file="${temp_prefix}_anno.tsv"         # Temp InterPro annotation
    temp_fa="${contig}.temp.fa"                      # Temp FASTA file for contig

    # Extract domain annotations for genes in the current contig
    awk -F '\t' -v name="${mcscan_contig}_" '$1 ~ name {print $1,$5,$6}' OFS='\t' "$interpro_anno/$strain_name.tsv" |
        sed 's/_/\t/3' | sed 's/^[^\t]*\t//' | awk '!seen[$0]++' > "$temp_anno_file"

    # Extract GFF entries and corresponding contig sequence
    awk -F '\t' -v name="${contig}" '$1 == name {print $0}'  "$orign_gff" > "$temp_gff_file"
    echo "$contig" > tempcontig
    seqkit grep -f tempcontig "$genome" | seqkit seq -w 50 > "$temp_fa"

    # Convert GFF and FASTA to GenBank format
    ./gff_to_gb_v10.pl -org "$strain_name" -gff "$temp_gff_file" -fasta "$temp_fa" > "$GBK/$i.gbk"

    # Annotate domains in the GenBank file
    python annotate_gbk_domain.py -i "$GBK/$i.gbk" -d "$temp_anno_file" -o temp_anno_gbk
    mv temp_anno_gbk "$GBK/$i.gbk"

    # Clean up temporary files
    rm -f "$temp_fa" "$temp_fa.index" tempcontig "$temp_gff_file" "$temp_anno_file"

done < "$node_list"


