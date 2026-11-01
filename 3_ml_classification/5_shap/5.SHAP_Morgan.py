from sklearn.ensemble import RandomForestClassifier
import shap
import pickle

# 1. Load data (saved in the previous step)
with open("./path/to/intermediate/data/4.Feature_Selection/rfe_results_Morgan.pkl", "rb") as f:
    rfe_results = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/X_train_filtered_Morgan.pkl", "rb") as f:
    X_train = pickle.load(f)
    
with open("./path/to/intermediate/data//4.Feature_Selection/y_train_filtered_Morgan.pkl", "rb") as f:
    y_train = pickle.load(f)

with open("./path/to/intermediate/data//4.Feature_Selection/feature_names_filtered_Morgan.pkl", "rb") as f:
    feature_names_filtered = pickle.load(f)
    
# 2. Optimal number of features (e.g., 10)
n_features_opt = 3
selected_features = rfe_results[n_features_opt]
matched_indices = [i for i, name in enumerate(feature_names_filtered) if name in selected_features]

# 3. Feature extraction
X_selected = X_train[:, matched_indices] # assumes NumPy array + column indices
  
# 4. Model training
clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_selected, y_train)

# 5. SHAP analysis
explainer = shap.Explainer(clf, X_selected)
shap_values = explainer(X_selected)

# 6. Plot (SHAP scores for important MACCS keys)
feature_names = [f"Feature {i}" for i in selected_features]

# Extract SHAP values for class 1 (positive side in binary classification)
shap_array_cluster1 = shap_values.values[:, :, 0]  # SHAP values for cluster 1
shap_array_cluster2 = shap_values.values[:, :, 1]  # SHAP values for cluster 2

# Select molecules with highest total SHAP contribution
total_shap_1 = np.abs(shap_array_cluster1).sum(axis=1)
top_indices_1 = np.argsort(-total_shap_1)[:10]  # top 10 molecules by SHAP contribution

# Visualization
shap_values_cluster1 = shap.Explanation(
    values = shap_array_cluster1,                         # shape: (24, 3)
    base_values = shap_values.base_values[:, 1],       # shape: (24,)
    data = shap_values.data
)

shap_values_cluster2 = shap.Explanation(
    values = shap_array_cluster2,                         # shape: (24, 3)
    base_values = shap_values.base_values[:, 0],       # shape: (24,)
    data = shap_values.data,
    feature_names = feature_names
)

shap.plots.bar(shap_values_cluster1)
shap.plots.bar(shap_values_cluster2, max_display=3)

shap.summary_plot(shap_values, X_selected)
