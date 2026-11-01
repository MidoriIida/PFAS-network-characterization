#############################
#   5 fold corss validation #
#############################
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import fetch_california_housing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pickle
import json

# 1. Load data (saved in the previous step)
with open("./path/to/intermediate/data/4.Feature_Selection/rfe_results_Morgan.pkl", "rb") as f:
    rfe_results = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/X_train_filtered_Morgan.pkl", "rb") as f:
    X_train = pickle.load(f)
    
with open("./path/to/intermediate/data//4.Feature_Selection/y_train_filtered_Morgan.pkl", "rb") as f:
    y_train = pickle.load(f)
    
with open("./path/to/intermediate/data//4.Feature_Selection/feature_names_filtered_Morgan.pkl", "rb") as f:
    feature_names_filtered = pickle.load(f)
    
# 2. Candidate feature counts
feature_counts = [1, 3, 5, 10, 20, 30, 50]

# 3. Run RFE and store results
rfe_5fcv_results = {}

for n_features in feature_counts: # cross-validate for each n_features
  selected_features = rfe_results[n_features]  # assumes RFE results were pre-saved
  matched_indices = [i for i, name in enumerate(feature_names_filtered) if name in selected_features]
  
  # 4. Extract by column name (DataFrame) or index (NumPy array)
  X_selected = X_train[:, matched_indices] # assumes NumPy array + column indices
    
  # 5. Define model
  clf = RandomForestClassifier(n_estimators=100, random_state=42)
  
  # 6. 5-fold cross-validation
  cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
  scores = cross_val_score(clf, X_selected, y_train, cv=cv, scoring='accuracy')
  
  # Print results
  print(f"5-fold CV Accuracy for n_features={n_features}:")
  print(f"Accuracy per fold: {np.round(scores, 3)}")
  print(f"Mean: {np.mean(scores):.4f} +/- {np.std(scores):.4f}")
  
  # Save if needed
  rfe_5fcv_results[n_features] = scores

# 7. Mean and variance of accuracy
mean_accuracies = []
for n_features in feature_counts:
  mean_accuracies.append(np.mean(rfe_5fcv_results[n_features]))
std_accuracies = []
for n_features in feature_counts:
  std_accuracies.append(np.std(rfe_5fcv_results[n_features]))

# 8. Plot
plt.figure(figsize=(8, 5))
plt.errorbar(feature_counts, mean_accuracies, yerr=std_accuracies, fmt='-o', capsize=5)
plt.xlabel("Number of Selected Features")
plt.ylabel("Cross-Validated Accuracy")
plt.title("Accuracy vs. Number of Selected Features (5-fold CV)")
plt.grid(True)
plt.tight_layout()
plt.show()

# 9. One-way ANOVA
stats.f_oneway(pd.Series(rfe_5fcv_results[1]),pd.Series(rfe_5fcv_results[3]), pd.Series(rfe_5fcv_results[5]), pd.Series(rfe_5fcv_results[10]),pd.Series(rfe_5fcv_results[20]), pd.Series(rfe_5fcv_results[30]), pd.Series(rfe_5fcv_results[50]))

# Define output file name
with open("./path/to/intermediate/data/4.Feature_Selection/rfe_5fcv_results_Morgan.pkl", "wb") as f:
    pickle.dump(rfe_5fcv_results, f)
    
# Expand dict to long format and save as CSV (UTF-8)
records = []
for key, array in rfe_5fcv_results.items():
    for idx, val in enumerate(array):
        records.append({"key": key, "index": idx, "value": val})

df = pd.DataFrame(records)
df.to_csv("./path/to/intermediate/data/4.Feature_Selection/rfe_5fcv_results_Morgan.csv", index=False)
