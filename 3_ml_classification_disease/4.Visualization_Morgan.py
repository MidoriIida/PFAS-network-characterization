##################################################################
#  Highlight SHAP-important Morgan fingerprint bits on          #
#  disease-space cluster molecules                              #
#  Important bits: 1398 (F-C-F-C-F), 1917 (C=O)               #
#  → Supplementary Figures S7, S8                              #
##################################################################
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pandas as pd

# 1. Load data
df = pd.read_csv("./path/to/intermediate/data/cas_to_smiles_disease.csv")

df_cl1 = df[df['label'] == 1]
df_cl2 = df[df['label'] == 2]

mols_cl1 = [Chem.MolFromSmiles(smi) for smi in df_cl1['SMILES'].tolist()]
mols_cl2 = [Chem.MolFromSmiles(smi) for smi in df_cl2['SMILES'].tolist()]
cas_cl1  = df_cl1['CAS'].tolist()
cas_cl2  = df_cl2['CAS'].tolist()

# 2. SHAP-important bits
important_bits = [1398, 1917]  # bit 1398: F-C-F-C-F, bit 1917: C=O

def get_highlight_lists(mols):
    atom_lists = []
    bond_lists = []
    for mol in mols:
        bit_info = {}
        AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, bitInfo=bit_info)
        highlight_atoms = set()
        highlight_bonds = set()
        for bit in important_bits:
            if bit in bit_info:
                for atom_idx, radius in bit_info[bit]:
                    env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)
                    highlight_bonds.update(env)
                    for bidx in env:
                        bond = mol.GetBondWithIdx(bidx)
                        highlight_atoms.add(bond.GetBeginAtomIdx())
                        highlight_atoms.add(bond.GetEndAtomIdx())
                    highlight_atoms.add(atom_idx)
        atom_lists.append(list(highlight_atoms))
        bond_lists.append(list(highlight_bonds))
    return atom_lists, bond_lists

hl_atoms_cl1, hl_bonds_cl1 = get_highlight_lists(mols_cl1)
hl_atoms_cl2, hl_bonds_cl2 = get_highlight_lists(mols_cl2)

# 3. Draw and save
img_cl1 = Draw.MolsToGridImage(
    mols_cl1, molsPerRow=4, subImgSize=(500, 500),
    legends=cas_cl1, highlightAtomLists=hl_atoms_cl1
)
img_cl2 = Draw.MolsToGridImage(
    mols_cl2, molsPerRow=4, subImgSize=(500, 500),
    legends=cas_cl2, highlightAtomLists=hl_atoms_cl2
)

out_dir = "/path/to/your/local/directory"
img_cl1.save(f"{out_dir}/FigureS7_disease_molecules_cl1.png")
img_cl2.save(f"{out_dir}/FigureS8_disease_molecules_cl2.png")

print(f"Cluster 1 (high-toxicity): {len(mols_cl1)} molecules → FigureS7")
print(f"Cluster 2 (low-toxicity):  {len(mols_cl2)} molecules → FigureS8")
