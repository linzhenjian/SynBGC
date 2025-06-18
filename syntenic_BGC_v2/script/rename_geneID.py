#!/usr/bin/env python3
"""
Usage
-----
python script.py input.gbk rename_list.txt renamed_output.gbk
"""

import sys
from pathlib import Path


def load_rename_dict(rename_file: str) -> dict:
    """Read two-column mapping file into a dict {old: new}."""
    rename_dict = {}
    with open(rename_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                old, new = parts
                rename_dict[old] = new
    return rename_dict


def rename_genes_and_contigs(input_gbk: str, output_gbk: str, rename_dict: dict):
    """
    • Replace /gene="OLD" with /gene="NEW"  (from rename_dict)
    • Remove underscores in contig IDs for every record
      – Detect the contig ID on the LOCUS line
      – Replace that exact ID (no other text) everywhere in the record
    """
    with open(input_gbk) as infile, open(output_gbk, "w") as outfile:
        current_old_ctg = current_new_ctg = None

        for line in infile:
            # --- (1) Handle LOCUS line: pick up record ID and strip underscores ---
            if line.startswith("LOCUS"):
                # LOCUS lines look like:  LOCUS       NW_012345678  12345 bp ...
                tokens = line.split()
                if len(tokens) >= 2:
                    current_old_ctg = tokens[1]
                    current_new_ctg = current_old_ctg.replace("_", "")
                    # Replace only the first occurrence of the contig ID on this line
                    line = line.replace(current_old_ctg, current_new_ctg, 1)

            # --- (2) Replace the contig ID elsewhere in the same record ---
            if current_old_ctg and current_old_ctg in line:
                line = line.replace(current_old_ctg, current_new_ctg)

            # --- (3) Rename gene qualifiers ---
            if '/gene="' in line:
                for old, new in rename_dict.items():
                    target = f'/gene="{old}"'
                    if target in line:
                        line = line.replace(target, f'/gene="{new}"')
                        break  # only one replacement per line

            # --- (4) Write the modified (or unmodified) line ---
            outfile.write(line)


def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py input.gbk rename_list.txt renamed_output.gbk", file=sys.stderr)
        sys.exit(1)

    input_gbk, rename_file, output_gbk = sys.argv[1:]
    if not Path(input_gbk).is_file():
        sys.exit(f"Error: {input_gbk} not found.")
    if not Path(rename_file).is_file():
        sys.exit(f"Error: {rename_file} not found.")

    rename_dict = load_rename_dict(rename_file)
    rename_genes_and_contigs(input_gbk, output_gbk, rename_dict)


if __name__ == "__main__":
    main()

