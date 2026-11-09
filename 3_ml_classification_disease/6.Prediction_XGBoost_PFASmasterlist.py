##################################################################
#  Apply disease-space XGBoost classifier to all PFAS in       #
#  the EPA CompTox PFAS Master List                            #
#  Output: Supplementary_Data_9_Classification_results_disease.txt
##################################################################
import pandas as pd
import numpy as np
import pickle

# 1. Load final model
with open("./path/to/intermediate/data/Final_xgboost_model_disease.pkl", "rb") as f:
    final_model = pickle.load(f)

# 2. Load master list features
with open("./path/to/intermediate/data/Fingerprint_Morgan_2048_chem_PFASmasterlist_disease.pkl", "rb") as f:
    data = pickle.load(f)

X_combined = data['X_combined']
df         = data['df']

# 3. Predict
y_pred  = final_model.predict(X_combined)
y_proba = final_model.predict_proba(X_combined)[:, 1]

# 4. Save results
df_out = df.copy()
df_out['predicted_label'] = y_pred
df_out['probability']     = y_proba
# label 1 = high-toxicity cluster, label 0 = low-toxicity cluster
df_out['cluster'] = df_out['predicted_label'].map({1: 'Cluster1_HighToxicity',
                                                    0: 'Cluster2_LowToxicity'})

df_out.to_csv(
    "./path/to/intermediate/data/Supplementary_Data_9_Classification_results_disease.txt",
    sep='\t', index=False
)

print(f"Predicted {len(df_out)} PFAS compounds")
print(f"  High-toxicity (Cluster 1): {(y_pred == 1).sum()}")
print(f"  Low-toxicity  (Cluster 2): {(y_pred == 0).sum()}")
print("Saved: Supplementary_Data_9_Classification_results_disease.txt")
