##################################################################
#  XGBoost classifier on Morgan fingerprints + chemical         #
#  descriptors using disease-space cluster labels (k=2)         #
#  Best params: n_estimators=300, max_depth=5, lr=0.3           #
#  CV accuracy: 0.772 ± 0.071                                   #
#  Outputs: Final_xgboost_model_disease.pkl                     #
##################################################################
import pandas as pd
import numpy as np
import pickle
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load combined features
with open("./path/to/intermediate/data/Fingerprint_Morgan_2048_chem_disease.pkl", "rb") as f:
    data = pickle.load(f)

X_combined = data['X_combined']
y = data['y']
df = data['df']

# Disease-space labels: 1 = high-toxicity cluster, 2 = low-toxicity cluster
y_binary = y - 1  # convert to 0/1

# 2. Train/test split
indices = np.arange(len(X_combined))
X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
    X_combined, y_binary, indices,
    test_size=0.3, random_state=42, stratify=y_binary
)

# 3. Best hyperparameters (from GridSearchCV in Analysis_251010)
best_params = {
    'colsample_bytree': 1.0,
    'learning_rate': 0.3,
    'max_depth': 5,
    'min_child_weight': 1,
    'n_estimators': 300,
    'subsample': 0.8,
    'random_state': 42,
    'eval_metric': 'logloss',
}

# 4. Train final model
final_model = XGBClassifier(**best_params)
final_model.fit(X_train, y_train)

cv_scores = cross_val_score(final_model, X_combined, y_binary, cv=5)
print(f"5-fold CV accuracy: {cv_scores.mean():.3f} (+/- {cv_scores.std():.3f})")
print(f"Train accuracy: {final_model.score(X_train, y_train):.3f}")
print(f"Test accuracy:  {final_model.score(X_test, y_test):.3f}")

# 5. Save model
with open("./path/to/intermediate/data/Final_xgboost_model_disease.pkl", "wb") as f:
    pickle.dump(final_model, f)

# 6. Evaluation
y_pred  = final_model.predict(X_test)
y_proba = final_model.predict_proba(X_test)[:, 1]
auc_score = roc_auc_score(y_test, y_proba)

print("\nClassification report:")
print(classification_report(y_test, y_pred,
                            target_names=['Low-toxicity (cl2)', 'High-toxicity (cl1)']))
print(f"ROC-AUC: {auc_score:.3f}")

# 7. Confusion matrix
cm = confusion_matrix(y_test, y_pred)
plt.figure(figsize=(6, 5))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
            xticklabels=['Low-tox', 'High-tox'],
            yticklabels=['Low-tox', 'High-tox'])
plt.title('Confusion Matrix (Disease-space classifier)')
plt.ylabel('True Label')
plt.xlabel('Predicted Label')
plt.tight_layout()
plt.savefig("./path/to/intermediate/data/confusion_matrix_disease.png", dpi=300)
plt.close()

# 8. ROC curve
fpr, tpr, _ = roc_curve(y_test, y_proba)
plt.figure(figsize=(6, 5))
plt.plot(fpr, tpr, label=f'XGBoost (AUC = {auc_score:.3f})', linewidth=2)
plt.plot([0, 1], [0, 1], 'k--', linewidth=1)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve - Disease-space classifier')
plt.legend()
plt.tight_layout()
plt.savefig("./path/to/intermediate/data/roc_curve_disease.png", dpi=300)
plt.close()

# 9. Save predictions
df_test = df.iloc[idx_test].copy()
df_test['true_label']      = y_test
df_test['predicted_label'] = y_pred
df_test['probability']     = y_proba
df_test['correct']         = (y_test == y_pred)
df_test.to_csv("./path/to/intermediate/data/predictions_XGBoost_disease.csv", index=False)

print("Saved: Final_xgboost_model_disease.pkl, confusion_matrix, roc_curve, predictions CSV")
