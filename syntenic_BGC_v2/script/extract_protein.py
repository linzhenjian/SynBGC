#!/usr/bin/env python3
import sys
from Bio import SeqIO
import warnings
from Bio import BiopythonParserWarning

warnings.simplefilter("ignore", BiopythonParserWarning)

def extract_proteins(gbk_file, output_fasta):
    with open(output_fasta, "w") as fasta_out:
        for record in SeqIO.parse(gbk_file, "genbank"):
            contig_id = record.id

            # -- organism name: keep only first two words ---------------
            raw_org = record.annotations.get("organism", "Unknown organism")
            first_two = "_".join(raw_org.split()[:2])  # e.g. "Escherichia_coli"
            organism = first_two

            for feature in record.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    gene_id = feature.qualifiers.get(
                        "gene",
                        feature.qualifiers.get("locus_tag", ["unknown"])
                    )[0]
                    protein_seq = feature.qualifiers["translation"][0]
                    header = f">{organism}_{contig_id}_{gene_id}"
                    fasta_out.write(f"{header}\n{protein_seq}\n")

# -- CLI ----------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_protein.py input.gbk output.fasta")
        sys.exit(1)

    extract_proteins(sys.argv[1], sys.argv[2])

