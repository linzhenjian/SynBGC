import argparse

def load_mapping(mapping_file):
    rename_dict = {}
    with open(mapping_file, 'r') as mf:
        for line in mf:
            names = line.strip().split(',')
            if len(names) > 1:
                primary = names[0]
                for name in names:
                    rename_dict[name] = primary
    return rename_dict

def rename_isoforms(data_file, output_file, rename_dict):
    with open(data_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                parts[0] = rename_dict.get(parts[0], parts[0])
                parts[1] = rename_dict.get(parts[1], parts[1])
            outfile.write('\t'.join(parts) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Unify gene isoform names using a mapping.')
    parser.add_argument('-m', '--map', required=True, help='Mapping file (CSV format)')
    parser.add_argument('-i', '--input', required=True, help='Input data file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    args = parser.parse_args()

    mapping = load_mapping(args.map)
    rename_isoforms(args.input, args.output, mapping)

if __name__ == '__main__':
    main()

