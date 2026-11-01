#############################
#   5 fold corss validation #
#############################
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import fetch_california_housing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pickle
import json

# 1. データの読み込み（前ステップで保存したテキスト）
with open("./path/to/intermediate/data/4.Feature_Selection/rfe_results_Morgan.pkl", "rb") as f:
    rfe_results = pickle.load(f)

with open("./path/to/intermediate/data/4.Feature_Selection/X_train_filtered_Morgan.pkl", "rb") as f:
    X_train = pickle.load(f)
    
with open("./path/to/intermediate/data//4.Feature_Selection/y_train_filtered_Morgan.pkl", "rb") as f:
    y_train = pickle.load(f)
    
with open("./path/to/intermediate/data//4.Feature_Selection/feature_names_filtered_Morgan.pkl", "rb") as f:
    feature_names_filtered = pickle.load(f)
    
#2. 特徴量数の候補
feature_counts = [1, 3, 5, 10, 20, 30, 50]

#3. RFE実行と結果格納
rfe_5fcv_results = {}

for n_features in feature_counts: #n_featuresを使って交差検証
  selected_features = rfe_results[n_features]  # 事前にRFE結果が保存されている前提
  matched_indices = [i for i, name in enumerate(feature_names_filtered) if name in selected_features]
  
  #4. データがDataFrameなら列名で抽出、NumPy配列なら番号で抽出
  X_selected = X_train[:, matched_indices] # NumPy配列 + 列番号と仮定
    
  #5. モデル定義
  clf = RandomForestClassifier(n_estimators=100, random_state=42)
  
  #6. 5-fold 交差検証
  cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
  scores = cross_val_score(clf, X_selected, y_train, cv=cv, scoring='accuracy')
  
  # 結果表示
  print(f"✅ 特徴量数 {n_features} における5-fold CV Accuracy:")
  print(f"各foldの正解率: {np.round(scores, 3)}")
  print(f"平均: {np.mean(scores):.4f} ± {np.std(scores):.4f}")
  
  # 保存しておく（必要なら）
  rfe_5fcv_results[n_features] = scores

#7. 精度の平均と分散
mean_accuracies = []
for n_features in feature_counts:
  mean_accuracies.append(np.mean(rfe_5fcv_results[n_features]))
std_accuracies = []
for n_features in feature_counts:
  std_accuracies.append(np.std(rfe_5fcv_results[n_features]))

#8. プロット
plt.figure(figsize=(8, 5))
plt.errorbar(feature_counts, mean_accuracies, yerr=std_accuracies, fmt='-o', capsize=5)
plt.xlabel("Number of Selected Features")
plt.ylabel("Cross-Validated Accuracy")
plt.title("Accuracy vs. Number of Selected Features (5-fold CV)")
plt.grid(True)
plt.tight_layout()
plt.show()

#9. 一限配置分散分析
stats.f_oneway(pd.Series(rfe_5fcv_results[1]),pd.Series(rfe_5fcv_results[3]), pd.Series(rfe_5fcv_results[5]), pd.Series(rfe_5fcv_results[10]),pd.Series(rfe_5fcv_results[20]), pd.Series(rfe_5fcv_results[30]), pd.Series(rfe_5fcv_results[50]))

# 保存先ファイル名を決めておく
with open("./path/to/intermediate/data/4.Feature_Selection/rfe_5fcv_results_Morgan.pkl", "wb") as f:
    pickle.dump(rfe_5fcv_results, f)
    
# 辞書を長形式に展開してcsvで保存（UTF-8で出力）
records = []
for key, array in rfe_5fcv_results.items():
    for idx, val in enumerate(array):
        records.append({"key": key, "index": idx, "value": val})

df = pd.DataFrame(records)
df.to_csv("./path/to/intermediate/data/4.Feature_Selection/rfe_5fcv_results_Morgan.csv", index=False)
