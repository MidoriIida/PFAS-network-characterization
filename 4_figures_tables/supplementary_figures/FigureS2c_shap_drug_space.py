##################################################
#   Figure S2c                                   #
#   Drug-space classifier                        #
#   SHAP bar chart                               #
#   Top 3 RFE features:                         #
#     Feature  650 (acidic oxygen: carboxylic or sulfonic acid) #
#     Feature 1088 (aromatic ring)               #
#     Feature 1398 (perfluoroalkyl chain)        #
#   Mean SHAP value by cluster (cluster 1/2)     #
#   モデル: X_train で学習                        #
#   SHAP:   全データ (X_all, n=35) で計算         #
#     → Feature 1088 がテストセットに偶然0個      #
#       だったため、全データで計算することで       #
#       代表性を確保する                          #
#   ターミナル(shap_env)で実行                    #
##################################################
import shap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import pickle

# ============================================================
# 1. データ読み込み（train + test）
# ============================================================

base = "/path/to/your/local/directory/Human_PPI_PFAS"

with open(f"{base}/Analysis_250412/4.Feature_Selection/X_train_filtered_Morgan.pkl", "rb") as f:
    X_train = pickle.load(f)
with open(f"{base}/Analysis_250412/4.Feature_Selection/y_train_filtered_Morgan.pkl", "rb") as f:
    y_train = pickle.load(f)
with open(f"{base}/Analysis_250412/4.Feature_Selection/X_test_filtered_Morgan.pkl", "rb") as f:
    X_test = pickle.load(f)
with open(f"{base}/Analysis_250412/4.Feature_Selection/y_test_filtered_Morgan.pkl", "rb") as f:
    y_test = pickle.load(f)
with open(f"{base}/Analysis_250412/4.Feature_Selection/rfe_results_Morgan.pkl", "rb") as f:
    rfe_results = pickle.load(f)
with open(f"{base}/Analysis_250412/4.Feature_Selection/feature_names_filtered_Morgan.pkl", "rb") as f:
    feature_names_filtered = pickle.load(f)

# ============================================================
# 2. n=3 の RFE 選択特徴量（Feature 650, 1088, 1398）を抽出
# ============================================================

n_features = 3
selected_features = rfe_results[n_features]
matched_indices = [i for i, name in enumerate(feature_names_filtered)
                   if name in selected_features]

X_sel_train = X_train[:, matched_indices]
X_sel_test  = X_test[:,  matched_indices]

# 全データ（n=35）: SHAP のグローバル解釈に使用
X_all = np.vstack([X_sel_train, X_sel_test])
y_all = np.concatenate([y_train, y_test])

print(f"Selected features (n={n_features}): {selected_features}")
print(f"Train: {X_sel_train.shape[0]}分子, Test: {X_sel_test.shape[0]}分子, All: {X_all.shape[0]}分子")
print(f"Feature 1088 in Train: {sum(X_sel_train[:,1]==1)}, Test: {sum(X_sel_test[:,1]==1)}, All: {sum(X_all[:,1]==1)}")

# ============================================================
# 3. RandomForest モデル構築（X_train で学習）
#    SHAP は全データ（X_all）で計算
# ============================================================

rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_sel_train, y_train)

explainer   = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(X_all)   # 全35分子で計算

# ============================================================
# 4. クラスターごとの平均 SHAP 値を計算
#    shap_values[:, :, 0] = class index 0 = label 1 (cluster 1) への寄与
# ============================================================

shap_df          = pd.DataFrame(shap_values[:, :, 0], columns=selected_features)
shap_df['label'] = np.array(y_all)

mean_shap = shap_df.groupby('label').mean().T
mean_shap.columns = ['cluster 1', 'cluster 2']
mean_shap['diff']  = mean_shap['cluster 1'] - mean_shap['cluster 2']
mean_shap          = mean_shap.sort_values('diff', ascending=True)  # 上: cluster1優位

print("\nMean SHAP by cluster:")
print(mean_shap[['cluster 1', 'cluster 2']].round(4))

# ============================================================
# 5. Bar chart 描画
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
