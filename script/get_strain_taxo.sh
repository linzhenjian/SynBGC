#!/bin/bash
#SBATCH --account=schmidt-np
#SBATCH --partition=schmidt-np

export TERM=xterm-256color
export PATH="/uufs/chpc.utah.edu/common/home/schmidt-group3/software/miniconda/bin:$PATH"

# Header
echo -e "class\torder\tfamily\tgenus\tspecies" > strain_taxo

# Loop through input file
while read A B; do
    id=$(echo "$B" | sed 's/_/ /')    # Replace underscores with spaces
    genus=$(echo "$id" | awk '{print $1}')  # First word = genus

    # Fetch XML once
    xml=$(esearch -db taxonomy -query "$id" < /dev/null | efetch -format xml)

    class=$(echo "$xml" | xtract -pattern LineageEx -block Taxon -if Rank -equals class -element ScientificName)
    order=$(echo "$xml" | xtract -pattern LineageEx -block Taxon -if Rank -equals order -element ScientificName)
    family=$(echo "$xml" | xtract -pattern LineageEx -block Taxon -if Rank -equals family -element ScientificName)

    # Output line
    printf "%s\t%s\t%s\t%s\t%s\n" "$class" "$order" "$family" "$genus" "$B" >> strain_taxo
done < "$1"

