import pandas as pd
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

# 1. Load data (saved in the previous step)
df = pd.read_csv("./path/to/intermediate/data/cas_to_smiles.csv")

# Create generator (bit vector, radius=2, 2048 bits)
generator = GetMorganGenerator(radius=2, fpSize=2048)

# 2. Generate Morgan fingerprints
fingerprints = []
valid_idx = []
for i, (cas, smi) in enumerate(zip(df['CAS'], df['SMILES'])):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        fp = generator.GetFingerprint(mol)
        arr = np.zeros((2048,), dtype=int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fingerprints.append(arr)
        valid_idx.append(i)

# 3. Keep only rows with valid fingerprints
df_valid = df.iloc[valid_idx].copy().reset_index(drop=True)
X = np.vstack(fingerprints)
y = df_valid['label'].astype(int).values

# 4. Save as pickle (X, y, df together for downstream use)
with open("./path/to/intermediate/data/Fingerprint_Morgan_2048.pkl", "wb") as f:
    pickle.dump({"X": X, "y": y, "df": df_valid}, f)

# 5. Also save as CSV
fp_df = pd.DataFrame(X)
fp_df["label"] = y
fp_df.to_csv("./path/to/intermediate/data/Fingerprint_Morgan_2048.csv", index=False, header=False)

print(f"done! {len(df_valid)} molecules, {X.shape[1]} bits")
