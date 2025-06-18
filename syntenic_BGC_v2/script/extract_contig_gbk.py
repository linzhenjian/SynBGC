from Bio import SeqIO
import sys

def extract_contig(input_gbk, output_gbk, contig_id):
    """Extracts a single contig from a multi-contig GenBank file."""
    with open(input_gbk) as gbk_file:
        records = [record for record in SeqIO.parse(gbk_file, "genbank") if record.id == contig_id]

    if records:
        SeqIO.write(records, output_gbk, "genbank")
        print(f"Extracted contig '{contig_id}' saved to {output_gbk}")
    else:
        print(f"Error: Contig '{contig_id}' not found in {input_gbk}")

# Example usage:
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_contig.py input.gbk output.gbk contig_id")
        sys.exit(1)
    extract_contig(sys.argv[1], sys.argv[2], sys.argv[3])

