##################################################################
#  Compute Morgan fingerprints (radius=2, 2048-bit) and         #
#  23 toxicity-related chemical descriptors for 44 PFAS         #
#  Outputs:                                                     #
#    Fingerprint_Morgan_2048_disease.pkl  {X, y, df}            #
#    Fingerprint_Morgan_2048_chem_disease.pkl  {X_combined, y}  #
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

# 1. Load SMILES
df = pd.read_csv("./path/to/intermediate/data/cas_to_smiles_disease.csv")

# 2. Morgan fingerprints
df['morgan_fp'] = df['SMILES'].apply(smiles_to_morgan_fp)
df_valid = df[df['morgan_fp'].notna()].copy()
X_fp = np.vstack(df_valid['morgan_fp'])
y = df_valid['label'].values

with open("./path/to/intermediate/data/Fingerprint_Morgan_2048_disease.pkl", "wb") as f:
    pickle.dump({'X': X_fp, 'y': y, 'df': df_valid}, f)

print(f"Morgan fingerprints: {X_fp.shape}")

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

scaler = StandardScaler()
X_desc_scaled = scaler.fit_transform(X_desc)

X_combined = np.hstack([X_fp, X_desc_scaled])
print(f"Chemical descriptors: {X_desc_scaled.shape}")
print(f"Combined features: {X_combined.shape}")

with open("./path/to/intermediate/data/Fingerprint_Morgan_2048_chem_disease.pkl", "wb") as f:
    pickle.dump({'X_combined': X_combined, 'y': y, 'df': df_valid,
                 'scaler': scaler, 'desc_names': DESC_NAMES}, f)

print("Saved: Fingerprint_Morgan_2048_disease.pkl, Fingerprint_Morgan_2048_chem_disease.pkl")
