import csv

input_file_path = 'no_errors.csv'
output_file_path = 'no_errors_clean.csv'

with open(input_file_path, mode='r', newline='') as input_file:
    csv_reader = csv.reader(input_file)
    rows = [row[2:] for row in csv_reader]

with open(output_file_path, mode='w', newline='') as output_file:
    csv_writer = csv.writer(output_file)
    csv_writer.writerows(rows)

print('Done!')
