from sklearn.ensemble import RandomForestClassifier
import shap
import pickle

# 1. データの読み込み（前ステップで保存したテキスト）
with open("./Analysis_250412/4.Feature_Selection/rfe_results_Morgan.pkl", "rb") as f:
    rfe_results = pickle.load(f)

with open("./Analysis_250412/4.Feature_Selection/X_train_filtered_Morgan.pkl", "rb") as f:
    X_train = pickle.load(f)
    
with open("./Analysis_250412//4.Feature_Selection/y_train_filtered_Morgan.pkl", "rb") as f:
    y_train = pickle.load(f)

with open("./Analysis_250412//4.Feature_Selection/feature_names_filtered_Morgan.pkl", "rb") as f:
    feature_names_filtered = pickle.load(f)
    
# 2. 最適な特徴数（例：10）
n_features_opt = 3
selected_features = rfe_results[n_features_opt]
matched_indices = [i for i, name in enumerate(feature_names_filtered) if name in selected_features]

# 3. 特徴量抽出
X_selected = X_train[:, matched_indices] # NumPy配列 + 列番号と仮定
  
# 4. モデル学習
clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_selected, y_train)

# 5. SHAP解析
explainer = shap.Explainer(clf, X_selected)
shap_values = explainer(X_selected)

# 6. プロット（重要なMACCSキーのSHAPスコア）
feature_names = [f"Feature {i}" for i in selected_features]

# クラス1（2クラス分類でpositive側）に対応する SHAP値だけを取り出す
shap_array_cluster1 = shap_values.values[:, :, 0]  # クラスター 1 に対するSHAP値
shap_array_cluster2 = shap_values.values[:, :, 1]  # クラスター 2 に対するSHAP値

# SHAPの総和が最大の分子を選ぶ
total_shap_1 = np.abs(shap_array_cluster1).sum(axis=1)
top_indices_1 = np.argsort(-total_shap_1)[:10]  # SHAP影響が大きい上位10個

# 可視化
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
