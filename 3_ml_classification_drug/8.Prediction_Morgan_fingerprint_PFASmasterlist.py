###############################
# Apply the trained model         #
#   and predict classification    #
###############################
import shap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import pickle
import seaborn as sns

### 1. load data ##########################
with open("./path/to/intermediate/data/4.Feature_Selection/X_train_filtered_Morgan.pkl", "rb") as f:
    X_train = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/y_train_filtered_Morgan.pkl", "rb") as f:
    y_train = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/rfe_results_Morgan.pkl", "rb") as f:
    rfe_results = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/feature_names_filtered_Morgan.pkl", "rb") as f:
    feature_names_filtered = pickle.load(f)

df = pd.read_csv("./path/to/intermediate/data/Fingerprint_Morgan_2048_pfas_masterlist.csv", header=None)

### 2. Prepare model ######################
n_features = 3 # number of features to select
selected_features = rfe_results[n_features]  # assumes RFE results were pre-saved
matched_indices = [i for i, name in enumerate(feature_names_filtered) if name in selected_features]
X_selected = X_train[:, matched_indices] # assumes NumPy array + column indices

### 3. Extract selected features from masterlist fingerprints #####
df_np = df.iloc[:, 0:2048].to_numpy()
df_selected = df_np[:, matched_indices]

### 4. Define model #######################
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_selected, y_train)  

### 5. Predict ###############################
predicted = rf.predict(df_selected) 
probabilities = rf.predict_proba(df_selected)

### 6. Organize and export results ###############
df_pred = df.iloc[:,2048:2051]
df_pred['prediction_label'] = pd.DataFrame(predicted)

df_prob = pd.DataFrame(probabilities)
df_pred['prediction_prob_cls1'] = df_prob.iloc[:, 0]
df_pred['prediction_prob_cls2'] = df_prob.iloc[:, 1]


df_feature = pd.DataFrame(df_selected)
df_pred = pd.concat([df_pred, df_feature], axis=1)

# Save full output (no header, for downstream R scripts)
output_path = "./path/to/intermediate/data/Prediction_pfas_masterlist.csv"
df_pred.to_csv(output_path, index=False, header=False)

# Save clean summary for merging (Supplementary Data 9)
df_summary = pd.DataFrame({
    'CASRN':             df.iloc[:, 2050].values,
    'SMILES':            df.iloc[:, 2049].values,
    'pred_drug_cluster': predicted,                    # 1 or 2
    'prob_drug_cl1':     df_prob.iloc[:, 0].values,
})
df_summary.to_csv(
    "./path/to/intermediate/data/Prediction_pfas_masterlist_drug.csv",
    index=False
)


