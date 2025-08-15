import argparse
from collections import defaultdict
from Bio import SeqIO

def load_domain_annotations(domain_file):
    domain_data = defaultdict(list)
    with open(domain_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                prot_id, domain_id, description = parts
                domain_data[prot_id].append(f"{domain_id}: {description}")
            elif len(parts) == 2:
                prot_id, domain_id = parts
                domain_data[prot_id].append(domain_id)
    return domain_data

def annotate_genbank(input_file, domain_data, output_file):
    records = []
    for record in SeqIO.parse(input_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                # Try multiple qualifiers to match the domain file
                pid = None
                if "protein_id" in feature.qualifiers:
                    pid = feature.qualifiers["protein_id"][0]
                elif "gene" in feature.qualifiers:
                    pid = feature.qualifiers["gene"][0]
                elif "locus_tag" in feature.qualifiers:
                    pid = feature.qualifiers["locus_tag"][0]

                if pid and pid in domain_data:
                    for note in domain_data[pid]:
                        feature.qualifiers.setdefault("note", []).append(f"domain: {note}")
        records.append(record)

    with open(output_file, "w") as out_handle:
        SeqIO.write(records, out_handle, "genbank")

def main():
    parser = argparse.ArgumentParser(description="Annotate GenBank file with conserved domain info.")
    parser.add_argument("-i", "--input", required=True, help="Input GenBank file")
    parser.add_argument("-d", "--domains", required=True, help="Tab-separated domain file")
    parser.add_argument("-o", "--output", required=True, help="Output annotated GenBank file")
    args = parser.parse_args()

    domain_data = load_domain_annotations(args.domains)
    annotate_genbank(args.input, domain_data, args.output)

if __name__ == "__main__":
    main()

