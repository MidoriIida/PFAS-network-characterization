##########################################
#   Figure S10                           #
#   Drug-space classifier                #
#   Boundary zone PFAS (prob 0.3-0.7)   #
#   Random sample of 20 molecules       #
#   Important bits: 650  (acidic oxygen) #
#                   1088 (aromatic ring) #
#                   1398 (F-C-F-C-F)    #
#   ターミナル(shap_env)で実行           #
##########################################
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pandas as pd
import random

random.seed(42)

# ============================================================
# 1. データ読み込み・境界域フィルタ
# ============================================================

df = pd.read_csv(
    "/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI_PFAS/Analysis_250412/Prediction_pfas_masterlist.csv",
    header=None
)
df.columns = ['url', 'SMILES', 'CASRN',
              'prediction_label', 'prob_cls1', 'prob_cls2',
              'c6', 'c7', 'c8']

# 境界域（予測確率 0.3–0.7）を抽出
boundary = df[(df['prob_cls1'] > 0.3) & (df['prob_cls1'] < 0.7)].copy()
boundary = boundary.dropna(subset=['SMILES'])

print(f"Boundary PFAS (0.3-0.7): {len(boundary)}")

# ============================================================
# 2. ランダムサンプリング（n=20）
# ============================================================

sample = boundary.sample(n=20, random_state=42).reset_index(drop=True)

cas_list = sample['CASRN'].tolist()
mols     = [Chem.MolFromSmiles(smi) for smi in sample['SMILES'].tolist()]

# SMILES が無効な場合は除外
valid = [(mol, cas) for mol, cas in zip(mols, cas_list) if mol is not None]
mols, cas_list = zip(*valid)
mols, cas_list = list(mols), list(cas_list)

print(f"Valid molecules: {len(mols)}")

# ============================================================
# 3. ハイライトするMorgan fingerprint ビット
#    bit  650: acidic oxygen (carboxylic or sulfonic acid)
#    bit 1088: aromatic ring
#    bit 1398: F-C-F-C-F（フルオロカーボン鎖）
# ============================================================

important_bits = [650, 1088, 1398]

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
# 4. グリッド描画・保存
# ============================================================

img = Draw.MolsToGridImage(
    mols,
    molsPerRow=4,
    subImgSize=(500, 500),
    legends=cas_list,
    highlightAtomLists=highlight_atom_lists
)

out_dir = "/Users/iida/Documents/論文/2024_Iida/第二稿/Code_for_Manuscript/Supplymentary_Figures"
img.save(f"{out_dir}/FigureS10_boundary_molecules_drug_space.png")
print(f"Saved: FigureS10_boundary_molecules_drug_space.png  ({len(mols)} molecules)")
