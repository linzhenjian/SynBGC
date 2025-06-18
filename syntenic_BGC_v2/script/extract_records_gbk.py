from Bio import SeqIO
import sys

def filter_gbk(input_gbk, output_gbk, gene_list_file):
    """Filters gene and CDS features from a GenBank file based on /gene qualifiers."""
    
    # Read list of gene names (from /gene qualifier)
    with open(gene_list_file) as f:
        gene_list = {line.strip() for line in f if line.strip()}

    records = []
    for record in SeqIO.parse(input_gbk, "genbank"):
        new_features = []

        for feature in record.features:
            if feature.type in {"gene", "CDS"}:
                # Only extract gene ID from the /gene qualifier
                gene_id = feature.qualifiers.get("gene", [None])[0]

                if gene_id and gene_id in gene_list:
                    new_features.append(feature)

        if new_features:
            record.features = new_features
            records.append(record)

    SeqIO.write(records, output_gbk, "genbank")
    print(f"Filtered GenBank file saved to {output_gbk}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_gbk.py input.gbk output.gbk gene_list.txt")
        sys.exit(1)
    filter_gbk(sys.argv[1], sys.argv[2], sys.argv[3])


