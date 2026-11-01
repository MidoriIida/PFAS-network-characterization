##########################
#  RFEで特徴量を刈り込み #
##########################
import pickle
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import numpy as np
import json

# 1. データの読み込み（前ステップで保存したテキスト）
with open("./path/to/intermediate/data/4.Feature_Selection/X_train_filtered_Morgan.pkl", "rb") as f:
    X = pickle.load(f)
    
with open("./path/to/intermediate/data//4.Feature_Selection/y_train_filtered_Morgan.pkl", "rb") as f:
    y = pickle.load(f)

with open("./path/to/intermediate/data//4.Feature_Selection/feature_names_filtered_Morgan.pkl", "rb") as f:
    feature_names = pickle.load(f)

# 2. 特徴量数の候補
feature_counts = [1, 3, 5, 10, 20, 30, 50]

# 3. 分類器（ランダムフォレストなど）
model = RandomForestClassifier(n_estimators=100, random_state=42) 

# 4. RFE実行と結果格納
rfe_results = {}

for n_features in feature_counts:
    if n_features > X.shape[1]:
        print(f"❌ スキップ: {n_features} は特徴量の数 {X.shape[1]} を超えています")
        continue

    selector = RFE(estimator=model, n_features_to_select=n_features, step=1)
    selector.fit(X, y)

    selected_mask = selector.support_
    selected_features = [f for f, keep in zip(feature_names, selected_mask) if keep]
    
    print(f"✅ 選ばれた特徴数: {n_features}")
    print("選ばれた特徴:", selected_features)
    print("-" * 40)

    # 保存しておく（必要なら）
    rfe_results[n_features] = selected_features

# 保存先ファイル名を決めておく
with open("./path/to/intermediate/data/4.Feature_Selection/rfe_results_Morgan.pkl", "wb") as f:
    pickle.dump(rfe_results, f)
    
# JSONファイルとして保存（UTF-8で出力）
with open("./path/to/intermediate/data/4.Feature_Selection/rfe_results_Morgan.json", "w", encoding="utf-8") as f:
    json.dump(rfe_results, f, ensure_ascii=False, indent=2)
    
