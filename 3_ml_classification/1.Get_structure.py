import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

# 1. Load data (saved in the previous step)
df = pd.read_csv("./path/to/intermediate/data/cas_to_smiles.csv")

# 2. Split data by cluster
df_cl1 = df[df['label'] == 1]
df_cl2 = df[df['label'] == 2]

# 3. Load multiple SMILES
smiles_list_sample_cl1 = df_cl1['SMILES'].dropna().tolist()
smiles_list_sample_cl2 = df_cl2['SMILES'].dropna().tolist()

# 4. Register to mol list
mols_cl1 = [Chem.MolFromSmiles(smi) for smi in smiles_list_sample_cl1]
mols_cl2 = [Chem.MolFromSmiles(smi) for smi in smiles_list_sample_cl2]

# 5. Draw structures in grid format
img_cl1 = Draw.MolsToGridImage(mols_cl1, molsPerRow=4, subImgSize=(250, 250), legends=df_cl1['CAS'].tolist())
img_cl2 = Draw.MolsToGridImage(mols_cl2, molsPerRow=4, subImgSize=(250, 250), legends=df_cl2['CAS'].tolist())

# 6. Save
img_path_cl1 = "./path/to/intermediate/data/smiles_structures_grid_cl1.png"
img_path_cl2 = "./path/to/intermediate/data/smiles_structures_grid_cl2.png"
img_cl1.save(img_path_cl1)
img_cl2.save(img_path_cl2)

