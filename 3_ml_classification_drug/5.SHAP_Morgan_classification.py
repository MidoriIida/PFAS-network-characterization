###############
# SHAP value calculation #
###############
import shap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import pickle
import seaborn as sns

### 1. Load data ##########################
with open("./path/to/intermediate/data/4.Feature_Selection/X_train_filtered_Morgan.pkl", "rb") as f:
    X_train = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/y_train_filtered_Morgan.pkl", "rb") as f:
    y_train = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/X_test_filtered_Morgan.pkl", "rb") as f:
    X_test = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/y_test_filtered_Morgan.pkl", "rb") as f:
    y_test = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/rfe_results_Morgan.pkl", "rb") as f:
    rfe_results = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/feature_names_filtered_Morgan.pkl", "rb") as f:
    feature_names_filtered = pickle.load(f)

### 2. Extract selected features ##########
n_features = 3
selected_features = rfe_results[n_features]
matched_indices = [i for i, name in enumerate(feature_names_filtered) if name in selected_features]
X_train_sel = X_train[:, matched_indices]
X_test_sel  = X_test[:,  matched_indices]

### 3. Train model on training data #######
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_train_sel, y_train)

### 4. Compute SHAP values on test data ###
explainer   = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(X_test_sel)  # list: [class0_array, class1_array]

print("Class 0 SHAP shape:", np.array(shap_values[0]).shape)
print("Class 1 SHAP shape:", np.array(shap_values[1]).shape)

### 5. Summary plots ######################
shap.summary_plot(shap_values[0], X_test_sel, feature_names=selected_features)  # cluster 1
shap.summary_plot(shap_values[1], X_test_sel, feature_names=selected_features)  # cluster 2

### 6. Mean SHAP values by cluster ########
shap_df = pd.DataFrame(shap_values[0], columns=selected_features)
shap_df['label'] = y_test

mean_shap_by_class = shap_df.groupby('label').mean().T
mean_shap_by_class.columns = ['cluster_1', 'cluster_2']
mean_shap_by_class['diff'] = mean_shap_by_class['cluster_1'] - mean_shap_by_class['cluster_2']
mean_shap_by_class = mean_shap_by_class.sort_values('diff', ascending=False)

mean_shap_by_class[['cluster_1', 'cluster_2']].plot(kind='barh', figsize=(6, 4))
plt.title('Average SHAP value (by cluster)')
plt.xlabel('SHAP value')
plt.tight_layout()
plt.show()

### 7. Scatter by cluster #################
fig, ax = plt.subplots(figsize=(10, 6))
for i, feature in enumerate(shap_df.columns[:-1]):
    for cls, color in zip([1, 2], ['blue', 'orange']):
        values = shap_df[shap_df['label'] == cls][feature]
        ax.scatter(values, [feature] * len(values), alpha=0.6,
                   label=cls if i == 0 else "", color=color)
ax.axvline(0, linestyle='--', color='gray')
ax.set_xlabel("SHAP value (impact on model output)")
ax.set_ylabel("Feature")
ax.set_title("SHAP summary plot by cluster")
ax.legend(title="Cluster")
plt.tight_layout()
plt.show()

### 8. Save ###############################
df_res = mean_shap_by_class[['cluster_1', 'cluster_2']]
df_res.to_csv("./path/to/intermediate/data/5.SHAP/SHAPvalue_Morgan_classification.csv", index=False)
