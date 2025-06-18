from Bio import SeqIO
import sys
import warnings
from Bio import BiopythonParserWarning

# Suppress BiopythonParserWarning
warnings.simplefilter("ignore", BiopythonParserWarning)

def gbk_to_gff3(gbk_file, output_file):
    with open(output_file, "w") as out:
        out.write("##gff-version 3\n")

        for record in SeqIO.parse(gbk_file, "genbank"):
            seqid = record.id
            for feature in record.features:
                # Always extract gene ID from the /gene qualifier.
                if feature.type == "gene":
                    gene_id = feature.qualifiers.get("gene", [f"gene_{feature.location.start}"])[0]
                    start = feature.location.start + 1  # convert to 1-based indexing for GFF3
                    end = feature.location.end
                    strand = '+' if feature.location.strand == 1 else '-' if feature.location.strand == -1 else '?'
                    out.write(f"{seqid}\tGenBank\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id}\n")

                elif feature.type == "mRNA":
                    # Use /gene attribute to extract the parent's gene id.
                    gene_id = feature.qualifiers.get("gene", [f"gene_{feature.location.start}"])[0]
                    # Generate transcript id from gene id (you may wish to modify this if there are multiple transcripts)
                    transcript_id = gene_id + ".mRNA"
                    start = feature.location.start + 1
                    end = feature.location.end
                    strand = '+' if feature.location.strand == 1 else '-' if feature.location.strand == -1 else '?'
                    out.write(f"{seqid}\tGenBank\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={transcript_id};Parent={gene_id}\n")

                elif feature.type == "CDS":
                    gene_id = feature.qualifiers.get("gene", [f"gene_{feature.location.start}"])[0]
                    transcript_id = gene_id + ".mRNA"
                    product = feature.qualifiers.get("product", [""])[0]
                    codon_start = int(feature.qualifiers.get("codon_start", [1])[0])
                    strand = '+' if feature.location.strand == 1 else '-' if feature.location.strand == -1 else '?'

                    # Handle possible compound locations
                    parts = feature.location.parts if hasattr(feature.location, "parts") else [feature.location]
                    phase = (3 - codon_start) % 3  # Initial phase for the first CDS part

                    for i, part in enumerate(parts):
                        start = part.start + 1
                        end = part.end
                        out.write(
                            f"{seqid}\tGenBank\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t"
                            f"ID={transcript_id}.CDS{i+1};Parent={transcript_id};product={product}\n"
                        )
                        phase = 0  # Only the first CDS part may have a non-zero phase

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gbk_to_gff3.py input.gbff output.gff3")
        sys.exit(1)
    else:
        gbk_to_gff3(sys.argv[1], sys.argv[2])

