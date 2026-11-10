##########################################
#   Figure S7                            #
#   Disease-space classifier             #
#   Molecular structure visualization   #
#   Cluster 1                            #
#   Important bits: 1398 (F-C-F-C-F)    #
#                   1917 (C=O)           #
#   Run in terminal (shap_env)            #
##########################################
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pandas as pd

# ============================================================
# 1. Load data (44 PFAS, disease-space classifier)
# ============================================================

df = pd.read_csv(
    "./path/to/intermediate/data/cas_to_smiles.csv"
)

df_cl1    = df[df['label'] == 1]
cas_list  = df_cl1['CAS'].dropna().tolist()
mols      = [Chem.MolFromSmiles(smi) for smi in df_cl1['SMILES'].dropna().tolist()]

print(f"Cluster 1: {len(mols)} molecules")

# ============================================================
# 2. Morgan fingerprint bits to highlight
#    bit 1398: F-C-F-C-F (fluorocarbon chain)
#    bit 1917: C=O (carbonyl group)
# ============================================================

important_bits = [1398, 1917]

highlight_atom_lists = []
for mol in mols:
    bit_info = {}
    AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=2048, bitInfo=bit_info
    )
    highlight_atoms = set()
    for bit in important_bits:
        if bit in bit_info:
            for atom_idx, radius in bit_info[bit]:
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)
                for bidx in env:
                    bond = mol.GetBondWithIdx(bidx)
                    highlight_atoms.add(bond.GetBeginAtomIdx())
                    highlight_atoms.add(bond.GetEndAtomIdx())
                highlight_atoms.add(atom_idx)
    highlight_atom_lists.append(list(highlight_atoms))

# ============================================================
# 3. Draw grid and save
# ============================================================

img = Draw.MolsToGridImage(
    mols,
    molsPerRow=4,
    subImgSize=(500, 500),
    legends=cas_list,
    highlightAtomLists=highlight_atom_lists
)

out_dir = "/path/to/your/local/directory/Supplymentary_Figures"
img.save(f"{out_dir}/FigureS7_disease_molecules_cl1.png")
print(f"Saved: FigureS7_disease_molecules_cl1.png")
