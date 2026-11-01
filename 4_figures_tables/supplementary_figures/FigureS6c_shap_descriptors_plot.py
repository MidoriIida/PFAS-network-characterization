#########################################
#　　モデルに寄与する重要な構造を特定   #
#########################################
#ターミナル(shap_env)で行う
import shap
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

# 1. データとモデルを読み込み
with open('/path/to/your/local/directory/Final_xgboost_model_Morgan_2048_chem.pkl', 'rb') as f:
    final_model = pickle.load(f)

with open('/path/to/your/local/directory/Fingerprint_Morgan_2048.pkl', 'rb') as f:
    data = pickle.load(f)

with open('/path/to/your/local/directory/Fingerprint_Morgan_2048_chem.pkl', 'rb') as f:
    data2 = pickle.load(f)

df = data['df']
X_combined = data2['X_combined']


# 2. SHAP Explainerを作成
explainer = shap.TreeExplainer(final_model)

# 3. SHAP値を計算（全データ）
shap_values = explainer.shap_values(X_combined)

# 4. SHAP値の形状を確認
print(f"SHAP値の形状: {shap_values.shape}") #SHAP値は (n_samples, n_features) の形状

# 5. 特徴量名のリストを作成
desc_names = [
    'MolWt', 'LogP', 'TPSA', 'HBD', 'HBA', 'RotBonds',
    'AromaticRings', 'AliphaticRings', 'NumRings', 'Heteroatoms',
    'ValenceElectrons', 'F', 'Cl', 'Br', 'I', 'O', 'N', 'S',
    'SpiroAtoms', 'BridgeheadAtoms', 'RadicalElectrons',
    'SaturatedRings', 'AromaticRings2'
] # 記述子の名前

n_fp = 2048
feature_names = [f'FP_{i}' for i in range(n_fp)] + desc_names

with open('./Analysis_251010/Final_xgboost_model_shap_Morgan_2048_chem.pkl', 'wb') as f:
  pickle.dump({'shap_values': shap_values, 'feature_names': feature_names}, f)
  
# ===== グローバルな重要度（平均絶対SHAP値） =====
print("\n" + "="*60)
print("SHAP値による特徴量重要度（グローバル）")
print("="*60)

mean_abs_shap = np.abs(shap_values).mean(axis=0)

# 記述子部分のSHAP値
desc_shap = mean_abs_shap[n_fp:]

print("\n記述子のSHAP重要度:")
desc_shap_sorted = sorted(zip(desc_names, desc_shap),
                          key=lambda x: x[1], reverse=True)
for i, (name, value) in enumerate(desc_shap_sorted, 1):
    print(f"{i:2d}. {name:20s}: {value:.6f}")

# Top 20特徴量
top_20_indices = np.argsort(mean_abs_shap)[-20:][::-1]
print("\n全特徴量のSHAP重要度 Top 20:")
for i, idx in enumerate(top_20_indices, 1):
    name = feature_names[idx]
    print(f"{i:2d}. {name:20s}: {mean_abs_shap[idx]:.6f}")

# ===== 記述子 Top 10 の SHAP Summary Plot =====
# 記述子部分のSHAP値（23次元）を抽出
shap_desc = shap_values[:, n_fp:]
X_desc    = X_combined[:, n_fp:]

# 平均絶対SHAP値で降順ソート → Top 10 インデックスを取得
mean_abs_desc = np.abs(shap_desc).mean(axis=0)
top10_idx     = np.argsort(mean_abs_desc)[::-1][:10]
top10_names   = [desc_names[i] for i in top10_idx]

print("\nTop 10 記述子（SHAP重要度順）:")
for rank, i in enumerate(top10_idx, 1):
    print(f"  {rank:2d}. {desc_names[i]:20s}: {mean_abs_desc[i]:.6f}")

# Top 10 に絞って summary_plot を出力
out_dir = "/path/to/your/local/directory/Supplymentary_Figures"

plt.figure(figsize=(8, 5))
shap.summary_plot(shap_desc[:, top10_idx], X_desc[:, top10_idx],
                  feature_names=top10_names, show=False)
plt.tight_layout()
out_path = f"{out_dir}/FigureS6c_shap_summary_plot.png"
plt.savefig(out_path, dpi=300, bbox_inches='tight')
print(f"\nSHAP Summary Plot (Top 10 descriptors) を保存: {out_path}")
plt.close()

# ===== Feature ImportanceとSHAPの比較 =====
print("\n" + "="*60)
print("Feature ImportanceとSHAP値の比較（記述子のみ）")
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

# 相関を計算
correlation = np.corrcoef(desc_importance, desc_shap)[0, 1]
print(f"\n相関係数: {correlation:.3f}")

# 比較プロット
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
plt.savefig('./Analysis_251010/importance_comparison.png', dpi=300, bbox_inches='tight')
print("\n比較プロットを保存: importance_comparison.png")
plt.close()

print("\n" + "="*60)
print("SHAP分析完了！")
print("="*60)

