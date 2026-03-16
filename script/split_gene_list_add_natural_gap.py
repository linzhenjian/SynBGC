import sys
import os
import re
import argparse
from determine_block_size import determine_block_size
from adjacent_numeric_gap import adjacent_numeric_gap
def read_gene_file(filename):
    """Reads a gene list file and returns a list of gene names."""
    with open(filename) as f:
        return [line.strip() for line in f if line.strip()]
def sum_insert_block_lengths(genes1, genes2):
    """
    Calculates the number of genes between consecutive anchors in genes1,
    based on their positions in genes2.
    """
    positions = {gene: i for i, gene in enumerate(genes2)}
    total_gap = 0

    for i in range(len(genes1) - 1):
        g1, g2 = genes1[i], genes1[i+1]
        if g1 in positions and g2 in positions:
            gap = positions[g2] - positions[g1] - 1
            total_gap += gap

    return total_gap

def find_large_insert_blocks(genes1, genes2, min_block_size):
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
    Extracts the numeric part from a gene ID like ..._g271.
    """
    match = re.search(r'_g(\d+)$', gene)
    return int(match.group(1)) if match else float('inf')

def find_large_gaps(genes, gap_threshold=26):
    """
    Finds positions where adjacent gene numbers differ by more than gap_threshold.
    """
    gap_positions = []
    for i in range(len(genes) - 1):
        num1 = extract_gene_number(genes[i])
        num2 = extract_gene_number(genes[i + 1])
        if abs(num2 - num1) > gap_threshold:
            gap_positions.append(i + 1)  # split after genes[i]
    return gap_positions

def split_gene_file(genes, split_positions, output_dir, output_base):
    """
    Splits gene list using given split positions.
    Always outputs at least one file, even if no splits occur.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    split_positions = sorted(split_positions)
    start = 0
    part_num = 1
    saved_files = []

    for end in split_positions + [len(genes)]:
        part = genes[start:end]
        # Always write the full file if no splits occurred
        if len(part) >= 3 or (start == 0 and end == len(genes)):
            outfile = os.path.join(output_dir, f"{output_base}_part{part_num}.split")
            with open(outfile, 'w') as f:
                f.write('\n'.join(part) + '\n')
            saved_files.append(outfile)
            part_num += 1
        start = end
    return saved_files

def main():
    """Main function with CLI arguments."""

    parser = argparse.ArgumentParser(description="Split gene lists based on insert blocks and gaps.")
    parser.add_argument("file1", help="First gene list file")
    parser.add_argument("file2", help="Second gene list file")
    parser.add_argument("outdir", help="Output directory")

    # New configurable thresholds
    parser.add_argument("--total_size", type=int, default=14,
                        help="Threshold for total insert block size (default: 14)")
    parser.add_argument("--min_median_gap", type=int, default=20,
                        help="Minimum median gap threshold (default: 20)")
    parser.add_argument("--max_median_gap", type=int, default=50,
                        help="Maximum median gap threshold (default: 50)")

    args = parser.parse_args()

    # DEFAULT_MEDIAN_GAP is identical to min_median_gap
    DEFAULT_MEDIAN_GAP = args.min_median_gap

    file1 = args.file1
    file2 = args.file2
    outdir = args.outdir.rstrip('/')
    outbase = os.path.splitext(os.path.basename(file1))[0]

    genes1 = read_gene_file(file1)
    genes2 = read_gene_file(file2)

    if len(genes1) <= 3:
        print("File1 has 3 or fewer lines. No segmentation will be done.")
        total_size = sum_insert_block_lengths(genes1, genes2)
        print(f"total_size is {total_size}")
        if total_size < args.total_size:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            outfile = os.path.join(outdir, f"{outbase}_part1.split")
            with open(outfile, 'w') as f:
                f.write('\n'.join(genes1) + '\n')
            print("Saved part: ", os.path.basename(outfile))
            return
        else:
            print("No file written because total_size >= threshold")
            return

    median_block_size = determine_block_size(file1, file2)
    print(f"Using median block size = {median_block_size} for splitting")
    insert_blocks = find_large_insert_blocks(genes1, genes2, min_block_size=median_block_size)
    print(f"Found {len(insert_blocks)} insert blocks >{median_block_size} genes")

    block_splits = find_split_positions(genes1, insert_blocks)
    print(f"Split positions from insert blocks: {block_splits}")

    if len(genes1) > 3:
        print(f"gene length = {len(genes1)}")
        median_gap = adjacent_numeric_gap(genes1)
        if median_gap is None:
            median_gap = DEFAULT_MEDIAN_GAP
            print(f"No natural gap found — using default median_gap = {DEFAULT_MEDIAN_GAP}")
        else:
            median_gap = median_gap - 1

        # Clamp median_gap between min and max
        if median_gap < args.min_median_gap:
            median_gap = args.min_median_gap
        if median_gap > args.max_median_gap:
            median_gap = args.max_median_gap

        print(f"Using median gap size = {median_gap} for splitting")
        gap_splits = find_large_gaps(genes1, gap_threshold=median_gap)
        print(f"Split positions from gene number gaps: {gap_splits}")
    else:
        gap_splits = []
        print("Skipped gap-based splitting (file1 has fewer than 3 genes).")

    split_positions = sorted(set(block_splits + gap_splits))
    print(f"Final split positions: {split_positions}")

    saved_files = split_gene_file(genes1, split_positions, outdir, outbase)

    if saved_files:
        print("Saved parts:")
        for f in saved_files:
            print("  ", os.path.basename(f))
    else:
        print("No output parts generated.")

if __name__ == "__main__":
    main()
