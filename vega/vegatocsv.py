import csv

# Read the .txt file
with open('report1.txt', 'r') as file:
    lines = file.readlines()

# Split lines by tabs and remove newline characters
data = [line.strip().split('\t') for line in lines]

# Write to a CSV file
with open('report1_table.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(data)