###################################
# Feature selection by correlation #
#     (remove redundant bits)      #
###################################
import numpy as np
import pickle
from sklearn.model_selection import train_test_split

# 1. Load fingerprint data
with open("./path/to/intermediate/data/Fingerprint_Morgan_2048.pkl", "rb") as f:
    data = pickle.load(f)

X = data["X"]
y = data["y"]
feature_names = np.array([f"Feature_{i}" for i in range(X.shape[1])])

# 2. Remove zero-variance features
std_series = X.std(axis=0)
mask = std_series > 0
X_filtered = X[:, mask]
feature_names_filtered = feature_names[mask]

print(f"features after zero-variance filter: {X_filtered.shape[1]} / {X.shape[1]}")

# 3. Remove highly correlated features (threshold = 0.7)
corr_matrix = np.corrcoef(X_filtered.T)
np.fill_diagonal(corr_matrix, 0)

n_features = X_filtered.shape[1]
selected = []
visited = set()

for i in range(n_features):
    if i in visited:
        continue
    group = set(np.where(corr_matrix[i] >= 0.7)[0])
    group.add(i)
    group_list = list(group)
    sub_corr = corr_matrix[np.ix_(group_list, group_list)]
    counts = (sub_corr >= 0.7).sum(axis=1)
    representative = group_list[np.argmax(counts)]
    selected.append(representative)
    visited.update(group)

selected_idx = sorted(set(selected))
X_reduced = X_filtered[:, selected_idx]
feature_names_selected = feature_names_filtered[selected_idx]

print(f"features after correlation filter: {X_reduced.shape[1]}")

# 4. Train / Test split (70/30)
indices = np.arange(X_reduced.shape[0])
X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
    X_reduced, y, indices, test_size=0.3, random_state=42
)

print(f"train: {X_train.shape[0]}, test: {X_test.shape[0]}")
print(f"test indices: {idx_test}")

# 5. Save as pickle
with open("./path/to/intermediate/data/4.Feature_Selection/X_train_filtered_Morgan.pkl", "wb") as f:
    pickle.dump(X_train, f)
with open("./path/to/intermediate/data/4.Feature_Selection/X_test_filtered_Morgan.pkl", "wb") as f:
    pickle.dump(X_test, f)
with open("./path/to/intermediate/data/4.Feature_Selection/y_train_filtered_Morgan.pkl", "wb") as f:
    pickle.dump(y_train, f)
with open("./path/to/intermediate/data/4.Feature_Selection/y_test_filtered_Morgan.pkl", "wb") as f:
    pickle.dump(y_test, f)
with open("./path/to/intermediate/data/4.Feature_Selection/feature_names_filtered_Morgan.pkl", "wb") as f:
    pickle.dump(feature_names_selected, f)

# 6. Also save as CSV (for downstream R scripts)
import pandas as pd
pd.DataFrame(X_train).to_csv(
    "./path/to/intermediate/data/4.Feature_Selection/X_train_filtered_Morgan.csv", index=False, header=False)
pd.DataFrame(X_test).to_csv(
    "./path/to/intermediate/data/4.Feature_Selection/X_test_filtered_Morgan.csv", index=False, header=False)
pd.DataFrame(y_train).to_csv(
    "./path/to/intermediate/data/4.Feature_Selection/Y_train_filtered_Morgan.csv", index=False, header=False)
pd.DataFrame(y_test).to_csv(
    "./path/to/intermediate/data/4.Feature_Selection/Y_test_filtered_Morgan.csv", index=False, header=False)
pd.DataFrame(feature_names_selected).to_csv(
    "./path/to/intermediate/data/4.Feature_Selection/feature_names_selected_Morgan.csv", index=False, header=False)

print("done!")
