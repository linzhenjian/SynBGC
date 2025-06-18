import sys

def read_gene_list(file_path):
    """Read genes from a file and return them as a list."""
    try:
        with open(file_path, "r") as file:
            genes = [line.strip() for line in file if line.strip()]
        return genes
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)

def filter_gene_list(genes):
    """Check if any two gene order numbers have a difference of 1."""
    gene_orders = sorted([int(g.split("_g")[-1]) for g in genes])  # Extract and sort gene order numbers
    
    # Check if any two consecutive numbers differ by exactly 1
    return any(gene_orders[i+1] - gene_orders[i] <= 2 for i in range(len(gene_orders) - 1))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py list1.txt")
        sys.exit(1)

    file_path = sys.argv[1]
    gene_list = read_gene_list(file_path)
    
    if filter_gene_list(gene_list):
        print(f"Valid file: {file_path}")
    else:
        print("No valid gene lists found.")

