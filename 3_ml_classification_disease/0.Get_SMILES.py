##############################################################
#  Get SMILES for 44 PFAS using disease-space cluster labels #
#  Reads list_cluster_disease.csv (k=2 from 2_disease_space) #
#  Output: cas_to_smiles_disease.csv                         #
##############################################################
import pandas as pd
import time
from pubchempy import get_compounds

def cas_to_cid(cas_number):
    try:
        compounds = get_compounds(cas_number, 'name')
        if compounds:
            return compounds[0].cid
        return None
    except Exception as e:
        print(f"  Error for {cas_number}: {e}")
        return None

def cid_to_smiles_pubchem(cid):
    try:
        compounds = get_compounds(cid, 'cid')
        if compounds:
            return compounds[0].isomeric_smiles
        return None
    except Exception as e:
        print(f"  Error for CID {cid}: {e}")
        return None

# 1. Load disease-space cluster labels
df = pd.read_csv("./path/to/intermediate/data/list_cluster_disease.csv")
df.reset_index(inplace=True)
df.columns = ['Unnamed: 0', 'CAS', 'label']

# 2. CAS → CID
cid_list = []
for i, cas in enumerate(df['CAS']):
    cid = cas_to_cid(str(cas).strip())
    cid_list.append(cid)
    print(f"[{i+1}/{len(df)}] CAS {cas} → CID {cid}")
    time.sleep(0.3)

df['CID'] = cid_list

# 3. CID → SMILES
smiles_list = []
for i, cid in enumerate(df['CID']):
    if pd.notna(cid):
        smi = cid_to_smiles_pubchem(int(cid))
        smiles_list.append(smi)
        print(f"[{i+1}/{len(df)}] CID {int(cid)} → SMILES obtained")
        time.sleep(0.3)
    else:
        smiles_list.append(None)

df['SMILES'] = smiles_list

# 4. Save
df_clean = df.dropna(subset=['SMILES'])
df_clean.to_csv("./path/to/intermediate/data/cas_to_smiles_disease.csv", index=False)

print(f"\nDone! {len(df_clean)}/{len(df)} SMILES obtained.")
print(f"Failed: {len(df) - len(df_clean)}")
