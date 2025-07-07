import argparse
from collections import defaultdict

def process_line(line):
    columns = line.strip().split()
    if not columns:
        return []
    family = columns[0]
    gene_entries = columns[1:]
    gene_dict = defaultdict(list)
    gene_order = []

    for gene_entry in gene_entries:
        parts = gene_entry.split('_')
        gene_name = '_'.join(parts[:3]) if len(parts) >= 3 else gene_entry
        if gene_name not in gene_dict:
            gene_order.append(gene_name)
        gene_dict[gene_name].append(gene_entry)

    output_lines = []
    for gene_name in gene_order:
        isoforms = ','.join(gene_dict[gene_name])
        output_lines.append([family, isoforms])

    return output_lines

def main():
    parser = argparse.ArgumentParser(description='Group gene isoforms together by gene name, separated by commas.')
    parser.add_argument('-i', '--input', required=True, help='Input file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    args = parser.parse_args()

    with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
        for line in infile:
            processed_lines = process_line(line)
            for entry in processed_lines:
                outfile.write('\t'.join(entry) + '\n')

if __name__ == '__main__':
    main()

