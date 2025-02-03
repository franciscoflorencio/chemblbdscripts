def split_file(input_file, lines_per_file=1410):
    with open(input_file, 'r') as file:
        file_number = 1
        lines = []
        for line in file:
            lines.append(line)
            if len(lines) == lines_per_file:
                with open(f'{input_file}_part{file_number}.txt', 'w') as output_file:
                    output_file.writelines(lines)
                lines = []
                file_number += 1
        if lines:
            with open(f'{input_file}_part{file_number}.txt', 'w') as output_file:
                output_file.writelines(lines)

if __name__ == "__main__":
    input_file = 'smiles_filtrado.txt'  # Replace with your file path
    split_file(input_file)