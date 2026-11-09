##################################################################
#  SHAP analysis for the disease-space XGBoost classifier       #
#  Run in shap_env (pip install shap)                           #
#  Outputs: shap_summary_disease.png, shap_descriptors_disease.png
##################################################################
import shap
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

# 1. Load model and features
with open("./path/to/intermediate/data/Final_xgboost_model_disease.pkl", "rb") as f:
    final_model = pickle.load(f)

with open("./path/to/intermediate/data/Fingerprint_Morgan_2048_chem_disease.pkl", "rb") as f:
    data = pickle.load(f)

X_combined = data['X_combined']
desc_names = data['desc_names']

# 2. Feature names
n_fp = 2048
feature_names = [f'FP_{i}' for i in range(n_fp)] + desc_names

# 3. SHAP values
explainer   = shap.TreeExplainer(final_model)
shap_values = explainer.shap_values(X_combined)
print(f"SHAP values shape: {shap_values.shape}")

with open("./path/to/intermediate/data/shap_disease.pkl", "wb") as f:
    pickle.dump({'shap_values': shap_values, 'feature_names': feature_names}, f)

# 4. Global importance (mean |SHAP|)
mean_abs_shap = np.abs(shap_values).mean(axis=0)

print("\nDescriptor SHAP importance:")
desc_shap = mean_abs_shap[n_fp:]
for name, val in sorted(zip(desc_names, desc_shap), key=lambda x: x[1], reverse=True):
    print(f"  {name:22s}: {val:.6f}")

print("\nTop 20 overall:")
for i, idx in enumerate(np.argsort(mean_abs_shap)[-20:][::-1], 1):
    print(f"  {i:2d}. {feature_names[idx]:22s}: {mean_abs_shap[idx]:.6f}")

# 5. SHAP summary plot (all features)
plt.figure(figsize=(10, 8))
shap.summary_plot(shap_values, X_combined, feature_names=feature_names,
                  max_display=20, show=False)
plt.tight_layout()
plt.savefig("./path/to/intermediate/data/shap_summary_disease.png", dpi=300, bbox_inches='tight')
plt.close()

# 6. SHAP plot for descriptors only
plt.figure(figsize=(10, 6))
shap.summary_plot(shap_values[:, n_fp:], X_combined[:, n_fp:],
                  feature_names=desc_names, show=False)
plt.tight_layout()
plt.savefig("./path/to/intermediate/data/shap_descriptors_disease.png", dpi=300, bbox_inches='tight')
plt.close()

print("Saved: shap_summary_disease.png, shap_descriptors_disease.png")
