from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import numpy as np

# 1. データ読み込み（前ステップで保存したもの）
df = pd.read_csv("./path/to/intermediate/data/cas_to_smiles.csv")

# SMILESとCAS列を取得
smiles_list = df['SMILES']
cas_list = df['CAS']

# ジェネレーターを作成（ビットベクトル、半径2、2048bit）
generator = GetMorganGenerator(radius=2, fpSize=2048)

# フィンガープリント生成
fingerprints = []
for cas, smi in zip(cas_list, smiles_list):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        fp = generator.GetFingerprint(mol)
        arr = np.zeros((2048,), dtype=int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fingerprints.append(arr)

# 6. ラベル取得（0 or 1）
y = df['label'].astype(int).values

# データフレームに変換
fp_df = pd.DataFrame(fingerprints)
fp_df['label'] = y  # 最後の列にクラスタラベルを追加

# 保存
output_path = "./path/to/intermediate/data/Fingerprint_Morgan_2048.csv"
fp_df.to_csv(output_path, index=False, header=False)
