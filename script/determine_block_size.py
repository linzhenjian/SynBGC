import re

def read_gene_file(filename):
    with open(filename) as f:
        return [line.strip() for line in f if line.strip()]

def find_large_insert_blocks(genes1, genes2, min_block_size=1):
    gene_set1 = set(genes1)
    insert_blocks = []
    current_block = []

    for gene in genes2:
        if gene not in gene_set1:
            current_block.append(gene)
        else:
            if len(current_block) > min_block_size:  # strictly greater
                insert_blocks.append(current_block[:])
            current_block = []
    if len(current_block) > min_block_size:
        insert_blocks.append(current_block)

    return insert_blocks

def extract_gene_number(gene):
    match = re.search(r'_g(\d+)$', gene)
    return int(match.group(1)) if match else float('inf')

def find_split_positions(genes1, insert_blocks):
    positions = []
    for block in insert_blocks:
        for gene in block:
            if gene in genes1:
                positions.append(genes1.index(gene))
                break
        else:
            for gene in genes1:
                if extract_gene_number(gene) > extract_gene_number(block[-1]):
                    positions.append(genes1.index(gene))
                    break
    return sorted(set(positions))

def split_gene_file(genes, split_positions):
    split_positions = sorted(split_positions)
    start = 0
    part_lengths = []

    for end in split_positions + [len(genes)]:
        part = genes[start:end]
        if len(part) > 0:
            part_lengths.append(len(part))
        start = end

    return part_lengths

def choose_median_block_size(results):
    """
    results: list of tuples (block_size, num_parts)
    """
    if not results:
        return None

    # Find max and min num_parts
    max_num_parts = max(results, key=lambda x: x[1])[1]
    min_num_parts = min(results, key=lambda x: x[1])[1]

    # First block size with min_num_parts
    first_min = next(b for b, n in results if n == min_num_parts)

    # Last block size with max_num_parts
    last_max = next(b for b, n in reversed(results) if n == max_num_parts)

    # Median logic
    candidates = sorted([first_min, last_max])
    if len(candidates) == 2:
        if candidates[1] - candidates[0] == 1:
            # If they are consecutive, pick the larger
            return candidates[1]
        else:
            # Otherwise take the true midpoint
            return (candidates[0] + candidates[1]) // 2
    return candidates[len(candidates)//2]

def determine_block_size(file1, file2, start_block=1, end_block=10, verbose=False):
    """
    Computes the median block size from file1 and file2.
    Returns an integer median block size.
    """
    genes1 = read_gene_file(file1)
    genes2 = read_gene_file(file2)
    results = []

    for block_size in range(start_block, end_block + 1):
        insert_blocks = find_large_insert_blocks(genes1, genes2, min_block_size=block_size)
        split_positions = find_split_positions(genes1, insert_blocks)
        part_lengths = split_gene_file(genes1, split_positions)
        num_parts = len(part_lengths)
        results.append((block_size, num_parts))

        if verbose:
            print(f"{block_size:<12} {num_parts:<12} {part_lengths}")

    median_block = choose_median_block_size(results)
    return median_block

# Allow running standalone
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python determine_block_size.py file1.txt file2.txt")
        sys.exit(1)

    median = determine_block_size(sys.argv[1], sys.argv[2], verbose=True)
    print("-" * 80)
    print(f"Median block size chosen: {median}")

