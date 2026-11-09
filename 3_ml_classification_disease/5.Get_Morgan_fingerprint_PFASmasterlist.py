##################################################################
#  Compute Morgan fingerprints + chemical descriptors for       #
#  all PFAS in the EPA CompTox PFAS Master List                 #
#  Input:  Chemical_List_pfasmaster.csv (from EPA CompTox)      #
#  Output: Fingerprint_Morgan_2048_chem_PFASmasterlist_disease.pkl
##################################################################
import pandas as pd
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Crippen, Lipinski
from sklearn.preprocessing import StandardScaler

def smiles_to_morgan_fp(smiles, radius=2, n_bits=2048):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            return np.array(fp)
        return None
    except:
        return None

def calculate_toxicity_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return [
        Descriptors.MolWt(mol),
        Crippen.MolLogP(mol),
        Descriptors.TPSA(mol),
        Lipinski.NumHDonors(mol),
        Lipinski.NumHAcceptors(mol),
        Lipinski.NumRotatableBonds(mol),
        Descriptors.NumAromaticRings(mol),
        Descriptors.NumAliphaticRings(mol),
        rdMolDescriptors.CalcNumRings(mol),
        Descriptors.NumHeteroatoms(mol),
        Descriptors.NumValenceElectrons(mol),
        smiles.count('F'),
        smiles.count('Cl'),
        smiles.count('Br'),
        smiles.count('I'),
        smiles.count('O'),
        smiles.count('N'),
        smiles.count('S'),
        rdMolDescriptors.CalcNumSpiroAtoms(mol),
        rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
        Descriptors.NumRadicalElectrons(mol),
        Lipinski.NumSaturatedRings(mol),
        Lipinski.NumAromaticRings(mol),
    ]

# 1. Load EPA PFAS Master List
# Download from: https://comptox.epa.gov/dashboard/chemical-lists/PFASMASTER
# Export → CSV, save and update path below
masterlist = pd.read_csv('/path/to/your/local/directory/Chemical_List_pfasmaster.csv')

# 2. Morgan fingerprints
masterlist['morgan_fp'] = masterlist['SMILES'].apply(
    lambda x: smiles_to_morgan_fp(x) if isinstance(x, str) else None
)
df_valid = masterlist[masterlist['morgan_fp'].notna()].copy()
X_fp = np.vstack(df_valid['morgan_fp'])
print(f"Master list: {len(masterlist)} total, {len(df_valid)} with valid SMILES")

# 3. Chemical descriptors
DESC_NAMES = [
    'MolWt', 'LogP', 'TPSA', 'HBD', 'HBA', 'RotBonds',
    'AromaticRings', 'AliphaticRings', 'NumRings', 'Heteroatoms',
    'ValenceElectrons', 'F', 'Cl', 'Br', 'I', 'O', 'N', 'S',
    'SpiroAtoms', 'BridgeheadAtoms', 'RadicalElectrons',
    'SaturatedRings', 'AromaticRings2',
]

descriptors = []
for smiles in df_valid['SMILES']:
    desc = calculate_toxicity_descriptors(smiles)
    descriptors.append(desc if desc is not None else [0] * len(DESC_NAMES))

X_desc = np.array(descriptors)

# 4. Load scaler fitted on training data
with open("./path/to/intermediate/data/Fingerprint_Morgan_2048_chem_disease.pkl", "rb") as f:
    train_data = pickle.load(f)
scaler = train_data['scaler']
X_desc_scaled = scaler.transform(X_desc)

X_combined = np.hstack([X_fp, X_desc_scaled])
print(f"Combined features: {X_combined.shape}")

# 5. Save
with open("./path/to/intermediate/data/Fingerprint_Morgan_2048_chem_PFASmasterlist_disease.pkl", "wb") as f:
    pickle.dump({'X_combined': X_combined, 'df': df_valid}, f)

print("Saved: Fingerprint_Morgan_2048_chem_PFASmasterlist_disease.pkl")
