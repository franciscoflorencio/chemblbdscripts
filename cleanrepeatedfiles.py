def filter_repeated_lines(input_file, output_file):
    seen_lines = set()
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line not in seen_lines:
                outfile.write(line)
                seen_lines.add(line)

filter_repeated_lines('smiles.txt', 'smiles_filtrado.txt')
        
