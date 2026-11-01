##########################################
#   Figure S8                            #
#   Disease-space classifier             #
#   Molecular structure visualization   #
#   Cluster 2                            #
#   Important bits: 1398 (F-C-F-C-F)    #
#                   1917 (C=O)           #
#   ターミナル(shap_env)で実行           #
##########################################
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pandas as pd

# ============================================================
# 1. データ読み込み（44 PFAS, disease-space classifier用）
# ============================================================

df = pd.read_csv(
    "/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI_PFAS/Analysis_251010/cas_to_smiles.csv"
)

df_cl2    = df[df['label'] == 2]
cas_list  = df_cl2['CAS'].dropna().tolist()
mols      = [Chem.MolFromSmiles(smi) for smi in df_cl2['SMILES'].dropna().tolist()]

print(f"Cluster 2: {len(mols)} molecules")

# ============================================================
# 2. ハイライトするMorgan fingerprint ビット
#    bit 1398: F-C-F-C-F（フルオロカーボン鎖）
#    bit 1917: C=O（カルボニル基）
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
# 3. グリッド描画・保存
# ============================================================

img = Draw.MolsToGridImage(
    mols,
    molsPerRow=4,
    subImgSize=(500, 500),
    legends=cas_list,
    highlightAtomLists=highlight_atom_lists
)

out_dir = "/Users/iida/Documents/論文/2024_Iida/第二稿/Code_for_Manuscript/Supplymentary_Figures"
img.save(f"{out_dir}/FigureS8_disease_molecules_cl2.png")
print(f"Saved: FigureS8_disease_molecules_cl2.png")
