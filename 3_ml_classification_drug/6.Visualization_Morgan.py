###########################
#  Visualize important structures  #
###########################
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pandas as pd

# 1. Load data (saved in the previous step)
df = pd.read_csv("./path/to/intermediate/data/cas_to_smiles.csv")

# 2. Split data by cluster
df_cl1 = df[df['label'] == 1]
df_cl2 = df[df['label'] == 2]
cas_list_1 = df_cl1['CAS'].dropna().tolist()
cas_list_2 = df_cl2['CAS'].dropna().tolist()

# 3. Load multiple SMILES
smiles_list_sample_cl1 = df_cl1['SMILES'].dropna().tolist()
smiles_list_sample_cl2 = df_cl2['SMILES'].dropna().tolist()

# 4. Register to mol list
mols_cl1 = [Chem.MolFromSmiles(smi) for smi in smiles_list_sample_cl1]
mols_cl2 = [Chem.MolFromSmiles(smi) for smi in smiles_list_sample_cl2]

# 5. Get important structural features
important_bits = [650, 1088, 1398]  # bit indices identified as important by SHAP

# 6. Generate drawing data
highlight_atom_lists_1 = []
highlight_atom_lists_2 = []
highlight_bond_lists_1 = []
highlight_bond_lists_2 = []

for mol in mols_cl1:
    bit_info = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, bitInfo=bit_info)

    highlight_atoms = set()
    highlight_bonds = set()
    for bit in important_bits:
        if bit in bit_info:
            for atom_idx, radius in bit_info[bit]:
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)
                highlight_bonds.update(env)
                atoms = set()
                for bidx in env:
                    bond = mol.GetBondWithIdx(bidx)
                    atoms.add(bond.GetBeginAtomIdx())
                    atoms.add(bond.GetEndAtomIdx())
                atoms.add(atom_idx)
                highlight_atoms.update(atoms)

    highlight_atom_lists_1.append(list(highlight_atoms))
    highlight_bond_lists_1.append(list(highlight_bonds))
    
for mol in mols_cl2:
    bit_info = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, bitInfo=bit_info)

    highlight_atoms = set()
    highlight_bonds = set()
    for bit in important_bits:
        if bit in bit_info:
            for atom_idx, radius in bit_info[bit]:
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)
                highlight_bonds.update(env)
                atoms = set()
                for bidx in env:
                    bond = mol.GetBondWithIdx(bidx)
                    atoms.add(bond.GetBeginAtomIdx())
                    atoms.add(bond.GetEndAtomIdx())
                atoms.add(atom_idx)
                highlight_atoms.update(atoms)

    highlight_atom_lists_2.append(list(highlight_atoms))
    highlight_bond_lists_2.append(list(highlight_bonds))

# 5. Draw structures in grid format
img_cl1 = Draw.MolsToGridImage(
    mols_cl1,
    molsPerRow=4,
    subImgSize=(500, 500),
    legends=cas_list_1,
    highlightAtomLists=highlight_atom_lists_1
)

img_cl2 = Draw.MolsToGridImage(
    mols_cl2,
    molsPerRow=4,
    subImgSize=(500, 500),
    legends=cas_list_2,
    highlightAtomLists=highlight_atom_lists_2
)

# 6. Save
img_path_cl1 = "./path/to/intermediate/data/highlighted_molecules_grid_cl1.png"
img_path_cl2 = "./path/to/intermediate/data/highlighted_molecules_grid_cl2.png"
img_cl1.save(img_path_cl1)
img_cl2.save(img_path_cl2)
