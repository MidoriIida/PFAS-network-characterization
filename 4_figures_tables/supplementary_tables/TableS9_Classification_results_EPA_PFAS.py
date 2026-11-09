##################################################################
#  Supplementary Data 9                                         #
#  Merge drug-space and disease-space predictions for all       #
#  PFAS in the EPA CompTox PFAS Master List                     #
#                                                               #
#  Prerequisites (run first):                                   #
#    3_ml_classification_drug/8.Prediction_Morgan_fingerprint_PFASmasterlist.py
#    3_ml_classification_disease/6.Prediction_XGBoost_PFASmasterlist.py
#                                                               #
#  Output: Supplementary_Data_9_Classification_results_for_all_EPA_PFAS.txt
##################################################################
import pandas as pd

# 1. Load drug-space predictions
# columns: CASRN, SMILES, pred_drug_cluster (1/2), prob_drug_cl1
# produced by: 3_ml_classification_drug/8.Prediction_Morgan_fingerprint_PFASmasterlist.py
df_drug = pd.read_csv("./path/to/intermediate/data/Prediction_pfas_masterlist_drug.csv")

# 2. Load disease-space predictions
# columns: CASRN, pred_disease_cluster (1/2), prob_disease_cl1
# produced by: 3_ml_classification_disease/6.Prediction_XGBoost_PFASmasterlist.py
df_disease = pd.read_csv("./path/to/intermediate/data/Prediction_pfas_masterlist_disease.csv")

# 3. Merge on CASRN
df = pd.merge(
    df_drug[['CASRN', 'SMILES', 'pred_drug_cluster', 'prob_drug_cl1']],
    df_disease[['CASRN', 'pred_disease_cluster', 'prob_disease_cl1']],
    on='CASRN', how='inner'
)

# 4. Assign group
def assign_group(row):
    d, s = row['pred_drug_cluster'], row['pred_disease_cluster']
    if   d == 1 and s == 1: return 'Group 1'
    elif d == 1 and s == 2: return 'Group 2'
    elif d == 2 and s == 1: return 'Group 3'
    else:                   return 'Group 4'

df['group'] = df.apply(assign_group, axis=1)

# 5. High-confidence flag
# TRUE when both classifiers agree with probability >= 0.7
# (prob_drug_cl1 >= 0.7 for drug=1, or <= 0.3 for drug=2; same for disease)
def is_high_confidence(row):
    drug_conf = row['prob_drug_cl1']    if row['pred_drug_cluster']    == 1 else 1 - row['prob_drug_cl1']
    dis_conf  = row['prob_disease_cl1'] if row['pred_disease_cluster'] == 1 else 1 - row['prob_disease_cl1']
    return drug_conf >= 0.7 and dis_conf >= 0.7

df['high_confidence'] = df.apply(is_high_confidence, axis=1)

# 6. Save
out_path = "./path/to/intermediate/data/Supplementary_Data_9_Classification_results_for_all_EPA_PFAS.txt"
df.to_csv(out_path, sep='\t', index=False)

print(f"Total PFAS predicted: {len(df)}")
print(df['group'].value_counts().sort_index())
print(f"\nHigh confidence: {df['high_confidence'].sum()} / {len(df)}")
print(f"Saved: {out_path}")
