from rdkit import Chem

# Function to canonicalize a SMILE
def canonicalize_smile(smile):
    mol = Chem.MolFromSmiles(smile)
    if mol:
        return Chem.MolToSmiles(mol, canonical=True)
    return None

# Read SMILES from input file
input_file = 'smiles.txt'
output_file = 'output_smiles.txt'

with open(input_file, 'r') as f:
    smiles = f.readlines()

# Strip whitespace and canonicalize each SMILE
canonical_smiles = set()
for smile in smiles:
    smile = smile.strip()
    canonical_smile = canonicalize_smile(smile)
    if canonical_smile:
        canonical_smiles.add(canonical_smile)

# Write unique canonical SMILES to output file
with open(output_file, 'w') as f:
    for smile in sorted(canonical_smiles):
        f.write(smile + '\n')

print(f"Canonical SMILES have been written to {output_file}")