#########################################
#  Identify structural features important to the model  #
#########################################
# Run in terminal (shap_env)
import shap
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

# 1. Load data and model
with open('./path/to/intermediate/data/Final_xgboost_model_Morgan_2048_chem.pkl', 'rb') as f:
    final_model = pickle.load(f)

with open('./path/to/intermediate/data/Fingerprint_Morgan_2048.pkl', 'rb') as f:
    data = pickle.load(f)

with open('./path/to/intermediate/data/Fingerprint_Morgan_2048_chem.pkl', 'rb') as f:
    data2 = pickle.load(f)

df = data['df']
X_combined = data2['X_combined']


# 2. Create SHAP Explainer
explainer = shap.TreeExplainer(final_model)

# 3. Compute SHAP values (all data)
shap_values = explainer.shap_values(X_combined)

# 4. Check SHAP value shape
print(f"SHAP value shape: {shap_values.shape}") # shape: (n_samples, n_features)

# 5. Create list of feature names
desc_names = [
    'MolWt', 'LogP', 'TPSA', 'HBD', 'HBA', 'RotBonds',
    'AromaticRings', 'AliphaticRings', 'NumRings', 'Heteroatoms',
    'ValenceElectrons', 'F', 'Cl', 'Br', 'I', 'O', 'N', 'S',
    'SpiroAtoms', 'BridgeheadAtoms', 'RadicalElectrons',
    'SaturatedRings', 'AromaticRings2'
] # descriptor names

n_fp = 2048
feature_names = [f'FP_{i}' for i in range(n_fp)] + desc_names

with open('./path/to/intermediate/data/Final_xgboost_model_shap_Morgan_2048_chem.pkl', 'wb') as f:
  pickle.dump({'shap_values': shap_values, 'feature_names': feature_names}, f)
  
# ===== Global importance (mean absolute SHAP value) =====
print("\n" + "="*60)
print("Feature importance by SHAP (global)")
print("="*60)

mean_abs_shap = np.abs(shap_values).mean(axis=0)

# SHAP values for descriptor features
desc_shap = mean_abs_shap[n_fp:]

print("\nDescriptor SHAP importance:")
desc_shap_sorted = sorted(zip(desc_names, desc_shap),
                          key=lambda x: x[1], reverse=True)
for i, (name, value) in enumerate(desc_shap_sorted, 1):
    print(f"{i:2d}. {name:20s}: {value:.6f}")

# Top 20 features
top_20_indices = np.argsort(mean_abs_shap)[-20:][::-1]
print("\nTop 20 SHAP feature importance (all features):")
for i, idx in enumerate(top_20_indices, 1):
    name = feature_names[idx]
    print(f"{i:2d}. {name:20s}: {mean_abs_shap[idx]:.6f}")

# ===== SHAP Summary Plot for Top 10 Descriptors =====
# Extract SHAP values for descriptor features (23 dimensions)
shap_desc = shap_values[:, n_fp:]
X_desc    = X_combined[:, n_fp:]

# Sort descending by mean absolute SHAP -> get Top 10 indices
mean_abs_desc = np.abs(shap_desc).mean(axis=0)
top10_idx     = np.argsort(mean_abs_desc)[::-1][:10]
top10_names   = [desc_names[i] for i in top10_idx]

print("\nTop 10 descriptors (by SHAP importance):")
for rank, i in enumerate(top10_idx, 1):
    print(f"  {rank:2d}. {desc_names[i]:20s}: {mean_abs_desc[i]:.6f}")

# Output summary_plot restricted to Top 10
out_dir = "/path/to/your/local/directory/Supplymentary_Figures"

plt.figure(figsize=(8, 5))
shap.summary_plot(shap_desc[:, top10_idx], X_desc[:, top10_idx],
                  feature_names=top10_names, show=False)
plt.tight_layout()
out_path = f"{out_dir}/FigureS6c_shap_summary_plot.png"
plt.savefig(out_path, dpi=300, bbox_inches='tight')
print(f"\nSHAP Summary Plot (Top 10 descriptors) saved: {out_path}")
plt.close()

# ===== Comparison of Feature Importance and SHAP =====
print("\n" + "="*60)
print("Feature Importance vs. SHAP comparison (descriptors only)")
print("="*60)

feature_importance = final_model.feature_importances_
desc_importance = feature_importance[n_fp:]

comparison_df = pd.DataFrame({
    'Feature': desc_names,
    'Feature_Importance': desc_importance,
    'SHAP_Importance': desc_shap
})
comparison_df = comparison_df.sort_values('SHAP_Importance', ascending=False)
print(comparison_df.to_string(index=False))

# Compute correlation
correlation = np.corrcoef(desc_importance, desc_shap)[0, 1]
print(f"\nCorrelation coefficient: {correlation:.3f}")

# Comparison plot
fig, ax = plt.subplots(1, 2, figsize=(14, 6))

# Feature Importance
comparison_df_sorted_fi = comparison_df.sort_values('Feature_Importance', ascending=True)
ax[0].barh(comparison_df_sorted_fi['Feature'], comparison_df_sorted_fi['Feature_Importance'])
ax[0].set_xlabel('Feature Importance')
ax[0].set_title('XGBoost Feature Importance')

# SHAP Importance
comparison_df_sorted_shap = comparison_df.sort_values('SHAP_Importance', ascending=True)
ax[1].barh(comparison_df_sorted_shap['Feature'], comparison_df_sorted_shap['SHAP_Importance'])
ax[1].set_xlabel('Mean |SHAP value|')
ax[1].set_title('SHAP Importance')

plt.tight_layout()
plt.savefig('./path/to/intermediate/data/importance_comparison.png', dpi=300, bbox_inches='tight')
print("\nComparison plot saved: importance_comparison.png")
plt.close()

print("\n" + "="*60)
print("SHAP analysis complete!")
print("="*60)

