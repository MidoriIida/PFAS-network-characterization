###############################
# 作ったモデルを使ってみたい  #
#   　 master listを          #
#  morganfingerprintへ変換    #
###############################
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

### 1. load data ##########################
masterlist = pd.read_csv('/path/to/your/local/directory/Chemical List pfasmaster-2025-05-08.csv')

### 2. smiles を Morgan fingerprintへ ########
# SMILESとCAS列を取得
smiles_list = masterlist['SMILES'].dropna().astype(str).tolist() #10986
cas_list    = masterlist['CASRN'].dropna().astype(str).tolist() #12039
DTXSID_list = masterlist['DTXSID'].dropna().astype(str).tolist() #12039

# ジェネレーターを作成（ビットベクトル、半径2、2048bit）
generator = GetMorganGenerator(radius=2, fpSize=2048)

# フィンガープリントベクトルを構築
fingerprints = [] #10974
valid_labels = [] #10974
valid_smiles = [] #10974
valid_cas = []

for smi, label, cas in zip(masterlist['SMILES'], masterlist['DTXSID'], masterlist['CASRN']):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fp = generator.GetFingerprint(mol)
        arr = np.zeros((2048,), dtype=int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fingerprints.append(arr)
        valid_labels.append(label)
        valid_smiles.append(smi)
        valid_cas.append(cas)
    except Exception as e:
        print(f"Skipped: {smi} → {e}")
        
# データフレームに変換
X = np.array(fingerprints)
y = np.array(valid_labels)

df_valid = pd.DataFrame(X)
df_valid['label'] = y
df_valid['SMILES'] = valid_smiles
df_valid['CASRN'] = valid_cas

# 保存
output_path = "./path/to/intermediate/data/Fingerprint_Morgan_2048_pfas_masterlist.csv"
df_valid.to_csv(output_path, index=False, header=False)

