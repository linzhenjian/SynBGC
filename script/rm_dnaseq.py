import argparse

def remove_dna_sequences_multi_contigs(input_file, output_file):
    writing = True

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            stripped = line.strip()

            if stripped == "//":
                outfile.write(line)
                writing = True  # Reset for next record
            elif stripped.startswith("ORIGIN"):
                outfile.write(line)
                writing = False  # Stop writing DNA sequence
            elif writing:
                outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove DNA sequences from multi-record GenBank file, keeping ORIGIN and // lines.")
    parser.add_argument("-i", "--input", required=True, help="Input GenBank (.gbk) file")
    parser.add_argument("-o", "--output", required=True, help="Output file")
    args = parser.parse_args()

    remove_dna_sequences_multi_contigs(args.input, args.output)

