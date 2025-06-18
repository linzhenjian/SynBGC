#!/usr/bin/env python3
import sys
from Bio import SeqIO

def extract_gene_locations(gbk_file, output_file):
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(gbk_file, "genbank"):
            contig = record.id

            # -- Extract only the first two words of the organism name --
            raw_org = record.annotations.get("organism", "Unknown organism")
            organism = "_".join(raw_org.split()[:2])  # e.g., "Escherichia_coli"

            org_contig = f"{organism}_{contig}"

            for feature in record.features:
                if feature.type == "gene":
                    gene_id = feature.qualifiers.get(
                        "gene",
                        feature.qualifiers.get("locus_tag", ["unknown"])
                    )[0]
                    full_id = f"{org_contig}_{gene_id}"

                    start = int(feature.location.start) + 1  # 1-based
                    end = int(feature.location.end)

                    out.write(f"{org_contig}\t{full_id}\t{start}\t{end}\n")

# -- CLI -----------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_gene_locations.py input.gbk output.tsv")
        sys.exit(1)

    extract_gene_locations(sys.argv[1], sys.argv[2])
