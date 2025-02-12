import pandas as pd

with_errors = 'with_errors.csv'

pd.read_csv(with_errors)['SMILES'].to_csv('smiles_errados.csv', index=False, header=False)