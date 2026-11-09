##################################################################
#  Apply disease-space XGBoost classifier to all PFAS in       #
#  the EPA CompTox PFAS Master List                            #
#  Output: Supplementary_Data_7_Classification_results_for_all_EPA_PFAS_disease.txt
##################################################################
import pandas as pd
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

# 4. Save results for merging (Supplementary Data 9)
# XGBoost labels: 1 = high-toxicity cluster (cl1), 0 = low-toxicity cluster (cl2)
# Convert to 1/2 to match drug-space convention
df_summary = pd.DataFrame({
    'CASRN':                df['CASRN'].values,
    'pred_disease_cluster': y_pred + 1,   # 0/1 → 1/2 (1=high-tox, 2=low-tox)
    'prob_disease_cl1':     y_proba,
})
df_summary.to_csv(
    "./path/to/intermediate/data/Prediction_pfas_masterlist_disease.csv",
    index=False
)

print(f"Predicted {len(df_summary)} PFAS compounds")
print(f"  High-toxicity (Cluster 1): {(y_pred == 1).sum()}")
print(f"  Low-toxicity  (Cluster 2): {(y_pred == 0).sum()}")
print("Saved: Prediction_pfas_masterlist_disease.csv")
