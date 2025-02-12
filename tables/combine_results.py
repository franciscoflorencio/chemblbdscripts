import os
import pandas as pd

input_dir = '../results'
output_csv = 'combined_results.csv'

combined_df = pd.DataFrame()

for filename in os.listdir(input_dir):
    if filename.endswith('.txt'):
        file_path = os.path.join(input_dir, filename)

        with open(file_path, 'r') as file:
            lines = file.readlines()

        for idx, line in enumerate(lines):
            if line.startswith("No."):
                data_start_idx = idx
                break

        header = lines[data_start_idx].strip().split('\t')
        data_lines = [line.strip().split('\t') for line in lines[data_start_idx + 1:] if line.strip()]

        df = pd.DataFrame(data_lines, columns=header)

        combined_df = pd.concat([combined_df, df], ignore_index=True)

combined_df.to_csv(output_csv, index=False)
print(f"Combined results saved to {output_csv}")