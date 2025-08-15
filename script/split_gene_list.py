import sys
import os

def read_gene_file(filename):
    """Reads a gene list file and returns a list of gene names."""
    with open(filename) as f:
        return [line.strip() for line in f if line.strip()]

def find_large_insert_blocks(genes1, genes2, min_block_size=3):
    """
    Identifies blocks of genes that are present in genes2 but missing from genes1.
    Only returns blocks longer than min_block_size.
    """
    gene_set1 = set(genes1)  # For fast lookup
    insert_blocks = []
    current_block = []

    for gene in genes2:
        if gene not in gene_set1:
            current_block.append(gene)  # Gene is new (not in file1)
        else:
            if len(current_block) > min_block_size:
                insert_blocks.append(current_block[:])  # Save block
            current_block = []  # Reset block
    if len(current_block) > min_block_size:
        insert_blocks.append(current_block)  # Catch trailing block

    return insert_blocks

def find_split_positions(genes1, insert_blocks):
    """
    For each insertion block, find the closest matching gene in file1,
    and use that index as a split point.
    """
    positions = []
    for block in insert_blocks:
        for gene in block:
            if gene in genes1:
                idx = genes1.index(gene)
                positions.append(idx)
                break
        else:
            for gene in genes1:
                if extract_gene_number(gene) > extract_gene_number(block[-1]):
                    idx = genes1.index(gene)
                    positions.append(idx)
                    break
    return sorted(set(positions))

def extract_gene_number(gene):
    """
    Helper function to extract the numeric part from a gene ID like ..._g271.
    Returns the gene number as an integer.
    """
    import re
    match = re.search(r'_g(\d+)$', gene)
    return int(match.group(1)) if match else float('inf')

def split_gene_file(genes, split_positions, output_dir, output_base):
    """
    Splits the input gene list using the given positions.
    Only writes parts with more than 3 genes.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    split_positions = sorted(split_positions)
    start = 0
    part_num = 1
    saved_files = []

    for end in split_positions + [len(genes)]:
        part = genes[start:end]
        if len(part) > 3 or (start == 0 and end == len(genes)):
            outfile = os.path.join(output_dir, f"{output_base}_part{part_num}.txt")
            with open(outfile, 'w') as f:
                f.write('\n'.join(part) + '\n')
            saved_files.append(outfile)
            part_num += 1
        start = end
    return saved_files

def main():
    """Main entry point: parses input arguments and runs processing."""
    if len(sys.argv) != 4:
        print("Usage: python split_gene_lists.py file1.txt file2.txt ./output/")
        return

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    outdir = sys.argv[3].rstrip('/')
    outbase = os.path.splitext(os.path.basename(file1))[0]

    genes1 = read_gene_file(file1)
    genes2 = read_gene_file(file2)

    if len(genes1) <= 10:
        print("File1 has 10 or fewer lines. No segmentation will be done.")
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outfile = os.path.join(outdir, f"{outbase}_part1.txt")
        with open(outfile, 'w') as f:
            f.write('\n'.join(genes1) + '\n')
        print("Saved part: ", os.path.basename(outfile))
        return

    insert_blocks = find_large_insert_blocks(genes1, genes2)
    print(f"Found {len(insert_blocks)} insert blocks >3 genes")

    split_positions = find_split_positions(genes1, insert_blocks)
    print(f"Splitting file1 at: {split_positions}")

    saved_files = split_gene_file(genes1, split_positions, outdir, outbase)

    if saved_files:
        print("Saved parts:")
        for f in saved_files:
            print("  ", os.path.basename(f))
    else:
        print("No output parts generated.")

if __name__ == "__main__":
    main()

