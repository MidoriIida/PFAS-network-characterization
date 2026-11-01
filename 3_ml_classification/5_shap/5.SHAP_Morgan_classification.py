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

### 1. load data ##########################
with open("./path/to/intermediate/data/4.Feature_Selection/X_test_filtered_Morgan.pkl", "rb") as f:
    X_test = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/y_test_filtered_Morgan.pkl", "rb") as f:
    y_test = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/rfe_results_Morgan.pkl", "rb") as f:
    rfe_results = pickle.load(f)

with open("./path/to/intermediate/data//4.Feature_Selection/feature_names_filtered_Morgan.pkl", "rb") as f:
    feature_names_filtered = pickle.load(f)

### 2. Extract data #######################
n_features = 3 # number of features to select
selected_features = rfe_results[n_features]  # assumes RFE results were pre-saved
matched_indices = [i for i, name in enumerate(feature_names_filtered) if name in selected_features]
X_selected = X_test[:, matched_indices] # assumes NumPy array + column indices

### 3. Define model #######################
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_selected, y_test)  
explainer = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(X_selected)
y_pred = rf.predict(X_selected)

print("Class 0 SHAP shape:", np.array(shap_values[0]).shape)
print("Class 1 SHAP shape:", np.array(shap_values[1]).shape)

# Compute SHAP values
shap_values = explainer.shap_values(X_selected)

# SHAP summary plot
shap_values_class1 = shap_values[:, :, 0] # SHAP values for class 1
shap.summary_plot(shap_values_class1, X_selected, feature_names=selected_features)

shap_values_class2 = shap_values[:, :, 1] # SHAP values for class 2 (positive)
shap.summary_plot(shap_values_class2, X_selected, feature_names=selected_features)

### 4. Visualize mean SHAP values ###################
# Convert class 1 SHAP values to DataFrame
shap_df = pd.DataFrame(shap_values[:, :, 0], columns=selected_features)

# Add labels
shap_df['label'] = y_test  # assumes y_test is a pandas Series or ndarray

# Compute mean SHAP values for class 0 and class 1
mean_shap_by_class = shap_df.groupby('label').mean().T
mean_shap_by_class.columns = ['cluster_1', 'cluster_2']
mean_shap_by_class['diff'] = mean_shap_by_class['cluster_1'] - mean_shap_by_class['cluster_2']
mean_shap_by_class = mean_shap_by_class.sort_values('diff', ascending=False)

mean_shap_by_class[['cluster_1', 'cluster_2']].plot(kind='barh', figsize=(6,4))
plt.title('Average SHAP value (class) ')
plt.xlabel('SHAP value')
plt.tight_layout()
plt.show()

# Overlay classification labels
fig, ax = plt.subplots(figsize=(10, 6))

for i, feature in enumerate(shap_df.columns[:-1]):  # last column is 'class'
    for cls, color in zip([1, 2], ['blue', 'orange']):
        values = shap_df[shap_df['label'] == cls][feature]
        ax.scatter(values, [feature] * len(values), alpha=0.6, label=cls if i == 0 else "", color=color)

ax.axvline(0, linestyle='--', color='gray')
ax.set_xlabel("SHAP value (impact on model output)")
ax.set_ylabel("Feature")
ax.set_title("SHAP summary plot by predicted class")
ax.legend(title="Predicted Class")
plt.tight_layout()
plt.show()

# Save
res = mean_shap_by_class[['cluster_1', 'cluster_2']]
df_res = pd.DataFrame(res)
df_res.to_csv("./path/to/intermediate/data/5.SHAP/SHAPvalue_Morgan_classification.csv", index=False)

