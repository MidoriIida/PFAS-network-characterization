from pubchempy import get_compounds
import pandas as pd
import time

# CAS → SMILES変換関数
def cas_to_smiles(cas_number):
    try:
        compound = get_compounds(cas_number, 'name')
        if compound:
            return compound[0].canonical_smiles
    except:
        return None

# 1. list_cluster.csv を読み込む（Rで出力したもの）
df = pd.read_csv("./Analysis_250412/list_cluster.csv")  # CAS, cluster
df.reset_index(inplace=True)
df.columns = ['Unnamed: 0','CAS', 'label']  # 列名を明示

# 3. SMILES列を追加（時間かかるので進捗確認しつつ）
smiles_list = []
for i, cas in enumerate(df['CAS']):
    smiles = cas_to_smiles(cas)
    smiles_list.append(smiles)
    print(f"[{i+1}/{len(df)}] {cas} → {smiles}")
    time.sleep(0.2)  # PubChemのアクセス制限対策として少し待つ

df['SMILES'] = smiles_list

# 4. SMILESが取得できた行だけにフィルター（オプション）
df_clean = df.dropna(subset=['SMILES'])

# 5. CSVに保存（変換結果）
df_clean.to_csv("./Analysis_250412/cas_to_smiles.csv", index=False)

print("done!")


