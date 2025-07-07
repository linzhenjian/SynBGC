#!/usr/bin/env python3
import argparse
import os  # Missing import
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description='Extract GFF entries for specific genes')
    parser.add_argument('-i', '--input', required=True, help='Input GFF file')
    parser.add_argument('-o', '--output', required=True, help='Output GFF file')
    parser.add_argument('-g', '--genes', required=True,
                      help='Comma-separated gene list (g1,g2,...) or file with one gene per line')
    return parser.parse_args()

def get_gene_ids(genes_arg):
    """Handle both comma-separated list and gene list files"""
    if os.path.isfile(genes_arg):
        with open(genes_arg, 'r') as f:
            return {line.strip() for line in f if line.strip()}
    # Fixed syntax below
    return {g.strip() for g in genes_arg.split(',')} if genes_arg else set()

def filter_gff():
    args = parse_arguments()
    target_genes = get_gene_ids(args.genes)
    
    with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            attributes = parts[8]
            gene_match = re.search(r'gene_id[= ]"?(?P<gene_id>\w+)"?', attributes)
            gene_id = gene_match.group('gene_id') if gene_match else None
            
            if not gene_id and parts[2] == 'gene':
                gene_id = attributes.strip()
            
            if gene_id in target_genes:
                outfile.write(line)

if __name__ == '__main__':
    filter_gff()
