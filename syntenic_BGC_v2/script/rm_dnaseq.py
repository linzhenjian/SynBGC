import argparse

def remove_dna_sequence_keep_origin_and_slash(input_file, output_file):
    writing = True
    origin_seen = False

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            stripped = line.strip()
            if stripped == "//":
                outfile.write(line)
                break
            elif stripped.startswith("ORIGIN"):
                origin_seen = True
                outfile.write(line)
                writing = False
            elif writing:
                outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove DNA sequence but keep ORIGIN and // lines from a GenBank file.")
    parser.add_argument("-i", "--input", required=True, help="Input GenBank (.gbk) file")
    parser.add_argument("-o", "--output", required=True, help="Output file")
    args = parser.parse_args()

    remove_dna_sequence_keep_origin_and_slash(args.input, args.output)

