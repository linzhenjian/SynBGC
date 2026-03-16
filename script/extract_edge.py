#!/usr/bin/env python3
import sys

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <species_file> <network_file> <out_file>")
    sys.exit(1)

species_file = sys.argv[1]
network_file = sys.argv[2]
out_file = sys.argv[3]

# Load species into a set
species = set()
with open(species_file) as f:
    for line in f:
        species.add(line.strip())

# Function to check if a value starts with any species in the set
def matches_species(value):
    for s in species:
        if value.startswith(s):
            return True
    return False

with open(network_file) as fin, open(out_file, "w") as fout:
    for line in fin:
        if line.startswith("Clustername"):  # Keep header
            fout.write(line)
            continue
        cols = line.strip().split("\t")
        if len(cols) < 2:
            continue
        col1, col2 = cols[0], cols[1]
        if matches_species(col1) and matches_species(col2):
            fout.write(line)

