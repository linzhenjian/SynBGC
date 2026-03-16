import sys
import os

def main():
    if len(sys.argv) != 4:
        print("Usage: script.py file1 file2 newfile")
        sys.exit(1)

    file1_path, file2_path, new_file_path = sys.argv[1:4]

    # Read input files
    with open(file1_path, 'r') as f:
        file1_genes = [line.strip() for line in f if line.strip()]
    if not file1_genes:
        print(f"Error: {file1_path} is empty")
        sys.exit(1)

    with open(file2_path, 'r') as f:
        file2_genes = [line.strip() for line in f if line.strip()]

    # Get first and last gene from file1
    first_gene, last_gene = file1_genes[0], file1_genes[-1]

    # Find positions in file2
    try:
        start_idx = file2_genes.index(first_gene)
        end_idx = file2_genes.index(last_gene)
    except ValueError as e:
        print(f"Error: {str(e)} not found in {file2_path}")
        sys.exit(1)

    # Get the section from file2
    subset_genes = file2_genes[start_idx:end_idx+1]

    # Add adjacent genes
    new_genes = []

    # Add genes from the section (restoring missing ones)
    new_genes.extend(subset_genes)

    # Ensure output directory exists (if directory is specified)
    output_dir = os.path.dirname(new_file_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Write output
    with open(new_file_path, 'w') as f:
        f.write('\n'.join(new_genes) + '\n')

    print(f"Processed file created at {new_file_path}")
    return 0

if __name__ == "__main__":
    sys.exit(main())

