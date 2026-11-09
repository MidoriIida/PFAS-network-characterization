##########################
#  Feature selection with RFE  #
##########################
import pickle
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import numpy as np
import json

# 1. Load data (saved in the previous step)
with open("./path/to/intermediate/data/4.Feature_Selection/X_train_filtered_Morgan.pkl", "rb") as f:
    X = pickle.load(f)
    
with open("./path/to/intermediate/data/4.Feature_Selection/y_train_filtered_Morgan.pkl", "rb") as f:
    y = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/feature_names_filtered_Morgan.pkl", "rb") as f:
    feature_names = pickle.load(f)

# 2. Candidate feature counts
feature_counts = [1, 3, 5, 10, 20, 30, 50]

# 3. Classifier (e.g., Random Forest)
model = RandomForestClassifier(n_estimators=100, random_state=42) 

# 4. Run RFE and store results
rfe_results = {}

for n_features in feature_counts:
    if n_features > X.shape[1]:
        print(f"Skipping: {n_features} exceeds the number of features {X.shape[1]}")
        continue

    selector = RFE(estimator=model, n_features_to_select=n_features, step=1)
    selector.fit(X, y)

    selected_mask = selector.support_
    selected_features = [f for f, keep in zip(feature_names, selected_mask) if keep]
    
    print(f"Selected feature count: {n_features}")
    print("Selected features:", selected_features)
    print("-" * 40)

    # Save if needed
    rfe_results[n_features] = selected_features

# Define output file name
with open("./path/to/intermediate/data/4.Feature_Selection/rfe_results_Morgan.pkl", "wb") as f:
    pickle.dump(rfe_results, f)
    
# Save as JSON (UTF-8 encoding)
with open("./path/to/intermediate/data/4.Feature_Selection/rfe_results_Morgan.json", "w", encoding="utf-8") as f:
    json.dump(rfe_results, f, ensure_ascii=False, indent=2)
    
