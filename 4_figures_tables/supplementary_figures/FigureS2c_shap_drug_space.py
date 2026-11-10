##################################################
#   Figure S2c                                   #
#   Drug-space classifier                        #
#   SHAP bar chart                               #
#   Top 3 RFE features:                         #
#     Feature  650 (acidic oxygen: carboxylic or sulfonic acid) #
#     Feature 1088 (aromatic ring)               #
#     Feature 1398 (perfluoroalkyl chain)        #
#   Mean SHAP value by cluster (cluster 1/2)     #
#   Model trained on X_train                      #
#   SHAP computed on all data (X_all, n=35)        #
#     Feature 1088 was absent in the test set,     #
#     so all data is used to ensure               #
#     representativeness.                          #
#   Run in terminal (shap_env)                    #
##################################################
import shap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import pickle

# ============================================================
# 1. Load data (train + test)
# ============================================================

base = "./path/to/intermediate/data/4.Feature_Selection"

with open(f"{base}/X_train_filtered_Morgan.pkl", "rb") as f:
    X_train = pickle.load(f)
with open(f"{base}/y_train_filtered_Morgan.pkl", "rb") as f:
    y_train = pickle.load(f)
with open(f"{base}/X_test_filtered_Morgan.pkl", "rb") as f:
    X_test = pickle.load(f)
with open(f"{base}/y_test_filtered_Morgan.pkl", "rb") as f:
    y_test = pickle.load(f)
with open(f"{base}/rfe_results_Morgan.pkl", "rb") as f:
    rfe_results = pickle.load(f)
with open(f"{base}/feature_names_filtered_Morgan.pkl", "rb") as f:
    feature_names_filtered = pickle.load(f)

# ============================================================
# 2. Extract RFE-selected features for n=3 (Feature 650, 1088, 1398)
# ============================================================

n_features = 3
selected_features = rfe_results[n_features]
matched_indices = [i for i, name in enumerate(feature_names_filtered)
                   if name in selected_features]

X_sel_train = X_train[:, matched_indices]
X_sel_test  = X_test[:,  matched_indices]

# All data (n=35): used for global SHAP interpretation
X_all = np.vstack([X_sel_train, X_sel_test])
y_all = np.concatenate([y_train, y_test])

print(f"Selected features (n={n_features}): {selected_features}")
print(f"Train: {X_sel_train.shape[0]} molecules, Test: {X_sel_test.shape[0]} molecules, All: {X_all.shape[0]} molecules")
print(f"Feature 1088 in Train: {sum(X_sel_train[:,1]==1)}, Test: {sum(X_sel_test[:,1]==1)}, All: {sum(X_all[:,1]==1)}")

# ============================================================
# 3. Build RandomForest model (trained on X_train)
#    SHAP computed on all data (X_all)
# ============================================================

rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_sel_train, y_train)

explainer   = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(X_all)   # computed on all 35 molecules

# ============================================================
# 4. Compute mean SHAP values per cluster
#    shap_values[:, :, 0] = class index 0 = contribution to label 1 (cluster 1)
# ============================================================

shap_df          = pd.DataFrame(shap_values[:, :, 0], columns=selected_features)
shap_df['label'] = np.array(y_all)

mean_shap = shap_df.groupby('label').mean().T
mean_shap.columns = ['cluster 1', 'cluster 2']
mean_shap['diff']  = mean_shap['cluster 1'] - mean_shap['cluster 2']
mean_shap          = mean_shap.sort_values('diff', ascending=True)  # top: cluster1 dominant

print("\nMean SHAP by cluster:")
print(mean_shap[['cluster 1', 'cluster 2']].round(4))

# ============================================================
# 5. Draw bar chart
# ============================================================

fig, ax = plt.subplots(figsize=(5, 3))

y_pos     = np.arange(len(mean_shap))
bar_width = 0.35

color_c1 = "#F8766D"   # salmon (cluster 1)
color_c2 = "#00BFC4"   # teal   (cluster 2)

bars1 = ax.barh(y_pos + bar_width / 2,
                mean_shap['cluster 1'],
                height=bar_width,
                color=color_c1, label='cluster 1')

bars2 = ax.barh(y_pos - bar_width / 2,
                mean_shap['cluster 2'],
                height=bar_width,
                color=color_c2, label='cluster 2')

ax.set_yticks(y_pos)
ax.set_yticklabels(mean_shap.index, fontsize=11)
ax.axvline(0, color='gray', linewidth=0.8, linestyle='--')
ax.set_xlabel("SHAP value", fontsize=11)
ax.legend(title="Class", fontsize=9, title_fontsize=9,
          loc='lower right')

plt.tight_layout()

out_dir = "/path/to/your/local/directory/Supplymentary_Figures"
plt.savefig(f"{out_dir}/FigureS2c_shap_drug_space.png", dpi=300, bbox_inches='tight')
print(f"\nSaved: FigureS2c_shap_drug_space.png")
plt.close()
