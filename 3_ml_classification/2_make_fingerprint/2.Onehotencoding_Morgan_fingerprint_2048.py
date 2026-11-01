from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import numpy as np

# 1. Load data (saved in the previous step)
df = pd.read_csv("./path/to/intermediate/data/cas_to_smiles.csv")

# Get SMILES and CAS columns
smiles_list = df['SMILES']
cas_list = df['CAS']

# Create generator (bit vector, radius=2, 2048 bits)
generator = GetMorganGenerator(radius=2, fpSize=2048)

# Generate fingerprints
fingerprints = []
for cas, smi in zip(cas_list, smiles_list):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        fp = generator.GetFingerprint(mol)
        arr = np.zeros((2048,), dtype=int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fingerprints.append(arr)

# 6. Get labels (0 or 1)
y = df['label'].astype(int).values

# Convert to data frame
fp_df = pd.DataFrame(fingerprints)
fp_df['label'] = y  # append cluster label as the last column

# Save
output_path = "./path/to/intermediate/data/Fingerprint_Morgan_2048.csv"
fp_df.to_csv(output_path, index=False, header=False)
