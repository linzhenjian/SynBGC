#! /usr/bin/python

import sys
import networkx as nx
from networkx import Graph

def load_edges(file_path):
    with open(file_path) as f:
        return [tuple(line.strip().split()) for line in f if line.strip()]

def get_connected_components(edge_list):
    g = Graph()
    g.add_edges_from(edge_list)
    components = list(nx.connected_components(g))
    # Sort by size, largest first
    return sorted(components, key=len, reverse=True)

def save_as_table(components, output_file):
    with open(output_file, 'w') as f:
        for i, comp in enumerate(components, start=1):
            gcf_id = f"GCF{i}"
            nodes = " ".join(sorted(comp))
            f.write(f"{gcf_id}\t{nodes}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_connected_components.py <input.txt> <output.tsv>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    edges = load_edges(input_path)
    components = get_connected_components(edges)
    save_as_table(components, output_path)

if __name__ == "__main__":
    main()


