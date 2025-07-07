import sys
import os

def parse_list_file(list_path):
    file_gene_map = {}
    with open(list_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 3:
                print(f"Skipping malformed line: {line}")
                continue
            filename, prefix, gene_orders = parts
            gene_names = [f"{prefix}_{gene.strip()}" for gene in gene_orders.split(",")]
            file_gene_map[filename] = gene_names
    return file_gene_map

def remove_gene_lines(file_gene_map):
    for filename, genes_to_remove in file_gene_map.items():
        if not os.path.isfile(filename):
            print(f"Warning: {filename} not found. Skipping.")
            continue

        with open(filename, "r") as infile:
            lines = infile.readlines()

        filtered_lines = [line for line in lines if line.strip() not in genes_to_remove]

        with open(filename, "w") as outfile:
            outfile.writelines(filtered_lines)

        print(f"Processed {filename}: removed {len(lines) - len(filtered_lines)} line(s).")

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py list")
        sys.exit(1)

    list_file = sys.argv[1]
    if not os.path.isfile(list_file):
        print(f"Error: List file '{list_file}' does not exist.")
        sys.exit(1)

    file_gene_map = parse_list_file(list_file)
    remove_gene_lines(file_gene_map)

if __name__ == "__main__":
    main()

