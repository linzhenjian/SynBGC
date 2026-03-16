#!/bin/bash
#SBATCH --account=schmidt-np
#SBATCH --partition=schmidt-np

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
    if [ -z "$class" ]; then
        class="unknown"
    fi
    if [ -z "$order" ]; then
        order="unknown"
    fi
    if [ -z "$family" ]; then
        family="unknown"
    fi
    # Output line
    printf "%s\t%s\t%s\t%s\t%s\n" "$class" "$order" "$family" "$genus" "$B" >> strain_taxo
done < "$1"

