##################################
# Build classification model using all features #
#        Morgan fingerprint       #
##################################
import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import classification_report, confusion_matrix

# 1. Load data (saved by 2.Onehotencoding_Morgan_fingerprint_2048.py)
with open("./path/to/intermediate/data/Fingerprint_Morgan_2048.pkl", "rb") as f:
    data = pickle.load(f)

X = data["X"]
y = data["y"]
df = data["df"]

print(f"n samples: {len(df)}")
print(f"cluster distribution:\n{df['label'].value_counts().sort_index()}")

# 2. Train/test split (stratified, 70/30)
indices = np.arange(len(X))
X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
    X, y, indices, test_size=0.3, random_state=42, stratify=y
)

print(f"\ntrain: {X_train.shape[0]}, test: {X_test.shape[0]}")

# 3. Train Random Forest
rf = RandomForestClassifier(
    n_estimators=200,
    max_depth=None,
    min_samples_split=2,
    min_samples_leaf=1,
    random_state=42,
    n_jobs=-1
)
rf.fit(X_train, y_train)

# 4. Evaluate
train_score = rf.score(X_train, y_train)
test_score  = rf.score(X_test,  y_test)
cv_scores   = cross_val_score(rf, X, y, cv=5)

print(f"\ntrain accuracy : {train_score:.3f}")
print(f"test accuracy  : {test_score:.3f}")
print(f"5-fold CV      : {cv_scores.mean():.3f} (+/- {cv_scores.std():.3f})")

y_pred = rf.predict(X_test)
print("\n=== Classification report ===")
print(classification_report(y_test, y_pred))
print("=== Confusion matrix ===")
print(confusion_matrix(y_test, y_pred))

# 5. Save model
with open("./path/to/intermediate/data/3.Constract_Full_Model/rf_model_Morgan.pkl", "wb") as f:
    pickle.dump(rf, f)
