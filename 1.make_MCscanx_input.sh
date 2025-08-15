#!/bin/bash
#SBATCH --account=schmidt-np    # account - abbreviated by -A
#SBATCH --partition=schmidt-np  # partition, abbreviated by -p

set -euo pipefail  # Exit on error, undefined variable, or pipeline failure

# ──────────────────────────────────────────────────────────────
# Help message displayed if no arguments or -h/--help is passed
# ──────────────────────────────────────────────────────────────
if [[ $# -eq 0 || "$1" == "-h" || "$1" == "--help" ]]; then
    cat <<EOF
synteny_BGC, Schmidt Lab, University of Utah

First, prepare a list of species names, Augustus GFF annotation files,
and genome assembly FASTA files, in this format:
Agelas_clathrodes    /path/to/augustus.gff    /path/to/genome.fna

Requirements: seqkit

Example:
  ./1.make_MCscanx_input.sh -a ./list_of_augustus_info -g ./folder_for_ncbi_gbk_files -o ./output/
EOF
    exit
fi

# ──────────────────────────────────────────────────────────────
# Get directory of current script, used for script paths later
# ──────────────────────────────────────────────────────────────
BIN_PATH="$(cd "$(dirname "$0")" && pwd)"
#BIN_PATH="$(dirname "$0")"
#BIN_PATH="/scratch/general/vast/zlin/sponge/REF_genome/syntenic_BGC"
# ──────────────────────────────────────────────────────────────
# Parse command-line options: -a (augustus list), -g (gbff dir), -o (output dir)
# ──────────────────────────────────────────────────────────────
while getopts "a:g:o:" opt; do
  case $opt in
    o) output="$OPTARG" ;;
    a) input_aug="$OPTARG" ;;
    g) input_gbk="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
  esac
done

# ──────────────────────────────────────────────────────────────
# At least one of -a or -g must be provided
# ──────────────────────────────────────────────────────────────
if [[ -z "${input_aug:-}" && -z "${input_gbk:-}" ]]; then
  echo "Error: Must provide -a (input_aug) or -g (input_gbk)" >&2
  exit 1
fi

# ──────────────────────────────────────────────────────────────
# Set output directory (default if not provided)
# ──────────────────────────────────────────────────────────────
if [[ -z "${output:-}" ]]; then
  output="./synteny_BGC_output"
fi
mkdir -p "$output"

# ──────────────────────────────────────────────────────────────
# Ensure required tools are available
# ──────────────────────────────────────────────────────────────
for cmd in awk seqkit sed mkdir; do
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "$cmd is required but not found. Aborting."
    exit 1
  }
done
#cores=$(getconf _NPROCESSORS_ONLN)
cores=$(nproc)
echo "paralle cores: $cores"
# ──────────────────────────────────────────────────────────────
# Create subdirectories for outputs
# ──────────────────────────────────────────────────────────────
mkdir -p "${output}/orign_prot" "${output}/orign_gff" "${output}/orign_gbk"

# ──────────────────────────────────────────────────────────────
# Generate gbk files if provided augutus files
# ──────────────────────────────────────────────────────────────
if [[ -f "${input_aug:-}" ]]; then
  # ──────────────────────────────────────────────────────────────
  # Function to convert GFF + genome FASTA to GenBank
  # Also extracts protein FASTA and cleaned GFF from GBK
  # ──────────────────────────────────────────────────────────────
  
  gff_2_gbk() {
    local organism="$1"
    local gff="$2"
    local genome="$3"
    #local gbk="$4"
    python "$BIN_PATH/script/gff_to_gbk.py" --organism "$organism" -g "$gff" -f "$genome" -o "$output/$organism.gbk"
    python "$BIN_PATH/script/remove_iso_multi_DNA.py" "$output/$organism.gbk" "$output/orign_gbk/${organism}.gbk"
    mv "$output/orign_gbk/${organism}.gbk"  "$output/$organism.gbk"
    python "$BIN_PATH/script/rm_dnaseq.py" -i "$output/${organism}.gbk" -o "$output/orign_gbk/${organism}.gbk"
    python "$BIN_PATH/script/extract_protein.py" "$output/orign_gbk/${organism}.gbk" "$output/orign_prot/${organism}.fa"
    sed -i 's/*//g' "$output/orign_prot/${organism}.fa"
    python "$BIN_PATH/script/extract_gff.py" "$output/orign_gbk/${organism}.gbk" "$output/orign_gff/${organism}.gff"
  }

  export -f gff_2_gbk
  export BIN_PATH output

  # Get number of CPU cores for parallel
  cores=$(getconf _NPROCESSORS_ONLN)
  echo "Using $cores cores"

 
  
  # Process each line of input_aug with GNU parallel
awk '
    NF >= 3 {
      organism = $1
      gff = $2
      genome = $3
      print organism, gff, genome
    }' "$input_aug" | awk '!seen[$1]++' \
| parallel -j "$cores" --colsep '[[:space:]]+' gff_2_gbk {1} {2} {3}
fi

# ──────────────────────────────────────────────────────────────
# Format NCBI .gbff files (if input_gbk directory is given)
# ──────────────────────────────────────────────────────────────

echo "[DEBUG] input_gbk = '${input_gbk:-}'"
if [[ -d "${input_gbk:-}" ]]; then
  echo "[DEBUG] input_gbk exists and is a directory"
else
  echo "[DEBUG] input_gbk is missing or not a directory"
fi


if [[ -d "${input_gbk:-}" ]]; then
  > $output/check_gbff
  find "$input_gbk" -maxdepth 1 -type f -name '*.gbff' -print0 | \
xargs -0 -I{} sh -c '
  gbff="{}"
  org=$(head -n 200 "$gbff" | awk '"'"'$1 == "ORGANISM" {print $2"_"$3}'"'"' | head -n 1)
  cds=$(head -n 4000 "$gbff" | awk '"'"'$1 == "CDS" {print $2}'"'"' | head -n 1)
  echo "$gbff $org $cds"
' >> "$output/check_gbff"

  awk '$3!="" {print $1,$2}' $output/check_gbff | awk '!seen[$2]++'  > $output/selected_gbff
  form_gbk () {
    local gbff="$1"
    local organism
    organism=$(awk '$1 == "ORGANISM" {print $2"_"$3}' "$gbff" | head -n 1)
    echo "$organism"

    python "$BIN_PATH/script/remove_iso_multi_DNA.py" "$gbff" "$output/${organism}.gbk"

    if [[ -s "$output/${organism}.gbk" ]]; then
      # File exists and is NOT empty — continue
      echo "[INFO] $output/${organism}.gbk is not empty, continuing..."

      python "$BIN_PATH/script/gbk_to_gff.py" "$output/${organism}.gbk" "$output/${organism}.gff"

      # Create gene name dictionary based on order
      awk -F '\t' '$3=="gene"' "$output/${organism}.gff" | \
        sed 's/>//g; s/<//g' | \
        sort -k1,1 -k4,4n | \
        nl | awk '{print $NF,"g"$1}' OFS='\t' | \
        sed 's/ID=//g' > "$output/${organism}.namedic"

      python "$BIN_PATH/script/rename_geneID.py" "$output/${organism}.gbk" "$output/${organism}.namedic" "$output/orign_gbk/${organism}.gbk"
      mv "$output/orign_gbk/${organism}.gbk"  "$output/${organism}.gbk"
      python "$BIN_PATH/script/rm_dnaseq.py" -i "$output/${organism}.gbk" -o "$output/orign_gbk/${organism}.gbk"

      python "$BIN_PATH/script/extract_protein.py" "$output/orign_gbk/${organism}.gbk" "$output/orign_prot/${organism}.fa"
      python "$BIN_PATH/script/extract_gff.py" "$output/orign_gbk/${organism}.gbk" "$output/orign_gff/${organism}.gff"

      # Clean up temp files
      rm "$output/${organism}.gff" "$output/${organism}.namedic" 
    else
      echo "[WARNING] $output/${organism}.gbk is empty. Skipping."
    fi
  }
  export -f form_gbk
  export BIN_PATH output

  # Run form_gbk in parallel over all *.gbff files
  awk '{print $1}' "$output/selected_gbff" | parallel -j "$cores" form_gbk
fi




#make comparison pairs from the gbk file names  

items=$(ls $output/orign_gbk/*.gbk | awk -F '/' '{print $NF}' | sed 's/[.]gbk//' | sort)
  set -- $items
  for i in $(seq 1 $#); do
    for j in $(seq $((i + 1)) $#); do
      eval "echo \${$i} \${$j}"
    done
  done > "$output/pairs"


