from pubchempy import get_compounds
import pandas as pd
import time

# CAS → SMILES conversion function
def cas_to_smiles(cas_number):
    try:
        compound = get_compounds(cas_number, 'name')
        if compound:
            return compound[0].canonical_smiles
    except:
        return None

# 1. Load list_cluster.csv (output from R)
df = pd.read_csv("./path/to/intermediate/data/list_cluster.csv")  # CAS, cluster
df.reset_index(inplace=True)
df.columns = ['Unnamed: 0','CAS', 'label']  # explicitly name columns

# 3. Add SMILES column (may take time; progress is printed)
smiles_list = []
for i, cas in enumerate(df['CAS']):
    smiles = cas_to_smiles(cas)
    smiles_list.append(smiles)
    print(f"[{i+1}/{len(df)}] {cas} → {smiles}")
    time.sleep(0.2)  # short delay to avoid PubChem rate limiting

df['SMILES'] = smiles_list

# 4. Keep only rows where SMILES was successfully retrieved (optional)
df_clean = df.dropna(subset=['SMILES'])

# 5. Save to CSV (conversion results)
df_clean.to_csv("./path/to/intermediate/data/cas_to_smiles.csv", index=False)

print("done!")


