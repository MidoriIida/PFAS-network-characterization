##########################
#   特徴量を削減してみる #
##########################
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.model_selection import cross_val_score
from xgboost import XGBClassifier

# 0. load data
with open('./Analysis_251010/Final_xgboost_model_shap_Morgan_2048_chem.pkl', 'rb') as f:
  data = pickle.load(f)

with open('./Analysis_251010/Fingerprint_Morgan_2048_chem.pkl', 'rb') as f:
  data2 = pickle.load(f)

shap_values = data['shap_values']
feature_names = data['feature_names']

X_combined = data2['X_combined']
y = data2['y']     # 1と2
y_binary = y - 1  # 0と1

best_params = {
    'colsample_bytree': 1.0,
    'learning_rate': 0.3,
    'max_depth': 5,
    'min_child_weight': 1,
    'n_estimators': 300,
    'subsample': 0.8,
    'random_state': 42,
    'eval_metric': 'logloss' #以下はwarningが出るため削除,'use_label_encoder': False
}

# 1. 特徴量数を変えた場合のCV正解率を計算
# 特徴量数のリストを定義
feature_counts = [5, 10, 20, 50, 100, 'all']

# shapの平均値を計算 
mean_abs_shap = np.abs(shap_values).mean(axis=0)

# 結果を保存する辞書
cv_results = {}


print("="*60)
print("特徴量数を変えた場合のCV正解率（各fold保存）")
print("="*60)

for n_features in feature_counts:
    if n_features == 'all':
        X_subset = X_combined
        n = X_combined.shape[1]
        label = f"All ({n})"
    else:
        top_indices = np.argsort(mean_abs_shap)[-n_features:]
        X_subset = X_combined[:, top_indices]
        n = n_features
        label = f"Top {n}"
    
    # CV実行
    clf_subset = XGBClassifier(**best_params)
    cv_scores = cross_val_score(clf_subset, X_subset, y_binary, cv=5)
    
    # 結果を保存
    cv_results[label] = cv_scores
    
    print(f"\n{label}:")
    print(f"  各fold: {[f'{s:.3f}' for s in cv_scores]}")
    print(f"  平均: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

# DataFrameに変換
cv_df = pd.DataFrame(cv_results)
print("\n" + "="*60)
print("CV正解率の一覧表")
print("="*60)
print(cv_df.round(3))

# CSVに保存
cv_df.to_csv('./Analysis_251010/cv_scores_by_features.csv', index=True)
print("\nCV結果を保存: cv_scores_by_features.csv")

# ===== t検定で統計的に比較 =====
print("\n" + "="*60)
print("t検定による統計的比較（vs 全特徴量）")
print("="*60)

baseline_label = list(cv_results.keys())[-1]  # 'All (2071)'
baseline_scores = cv_results[baseline_label]

comparison_results = []

for label, scores in cv_results.items():
    if label == baseline_label:
        continue
    
    # ペアt検定
    t_stat, p_value = stats.ttest_rel(baseline_scores, scores)
    
    mean_diff = baseline_scores.mean() - scores.mean()
    
    print(f"\n{label} vs {baseline_label}:")
    print(f"  平均差: {mean_diff:.4f}")
    print(f"  t統計量: {t_stat:.3f}")
    print(f"  p値: {p_value:.4f}")
    
    if p_value < 0.001:
        significance = "*** (p<0.001)"
    elif p_value < 0.01:
        significance = "** (p<0.01)"
    elif p_value < 0.05:
        significance = "* (p<0.05)"
    else:
        significance = "n.s. (p≥0.05)"
    
    print(f"  統計的有意性: {significance}")
    
    comparison_results.append({
        'Comparison': f"{label} vs {baseline_label}",
        'Mean_Diff': mean_diff,
        't_stat': t_stat,
        'p_value': p_value,
        'Significance': significance
    })

# 比較結果をDataFrameに
comparison_df = pd.DataFrame(comparison_results)
comparison_df.to_csv('./Analysis_251010/t_test_results.csv', index=False)
print("\nt検定結果を保存: t_test_results.csv")

# ===== ボックスプロット =====
print("\n" + "="*60)
print("ボックスプロットの作成")
print("="*60)

# データを整形（long format）
plot_data = []
for label, scores in cv_results.items():
    for score in scores:
        plot_data.append({'Features': label, 'CV Score': score})

plot_df = pd.DataFrame(plot_data)

# ボックスプロット
fig, ax = plt.subplots(figsize=(12, 6))

# seabornでボックスプロット
sns.boxplot(data=plot_df, x='Features', y='CV Score', ax=ax)
sns.swarmplot(data=plot_df, x='Features', y='CV Score', ax=ax, 
              color='black', alpha=0.5, size=6)

# 平均値をラインで追加
means = [cv_results[label].mean() for label in cv_results.keys()]
ax.plot(range(len(means)), means, 'r--', marker='D', 
        markersize=8, linewidth=2, label='Mean')

ax.set_ylabel('CV Accuracy', fontsize=12)
ax.set_xlabel('Number of Features', fontsize=12)
ax.set_title('Cross-Validation Accuracy by Number of Features', fontsize=14)
ax.legend()
ax.grid(axis='y', alpha=0.3)
ax.set_ylim(0.5, 1.0)

plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('./Analysis_251010/cv_boxplot.png', dpi=300, bbox_inches='tight')
print("ボックスプロットを保存: cv_boxplot.png")
plt.close()

# ===== バイオリンプロット（より詳細） =====
fig, ax = plt.subplots(figsize=(12, 6))

sns.violinplot(data=plot_df, x='Features', y='CV Score', ax=ax)
sns.swarmplot(data=plot_df, x='Features', y='CV Score', ax=ax, 
              color='black', alpha=0.6, size=6)

ax.set_ylabel('CV Accuracy', fontsize=12)
ax.set_xlabel('Number of Features', fontsize=12)
ax.set_title('Cross-Validation Accuracy Distribution by Number of Features', fontsize=14)
ax.grid(axis='y', alpha=0.3)
ax.set_ylim(0.5, 1.0)

plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('./Analysis_251010/cv_violinplot.png', dpi=300, bbox_inches='tight')
print("バイオリンプロットを保存: cv_violinplot.png")
plt.close()

# ===== 統計的サマリー =====
print("\n" + "="*60)
print("統計的サマリー")
print("="*60)

summary_stats = cv_df.describe().T
summary_stats['median'] = cv_df.median()
summary_stats = summary_stats[['mean', 'median', 'std', 'min', 'max']]

print(summary_stats.round(4))
summary_stats.to_csv('./Analysis_251010/cv_summary_stats.csv')
print("\nサマリー統計を保存: cv_summary_stats.csv")

print("\n分析完了！")
