##################################
# Build classification model using all features #
# 　 Morgan finger print         #
##################################
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
import shap
import matplotlib.pyplot as plt
import pickle

# 1. Load data (saved in the previous step)
with open("./path/to/intermediate/data/3.Constract_Full_Model//X_train_Morgan.pkl", "rb") as f:
    X_train = pickle.load(f)

with open("./path/to/intermediate/data/3.Constract_Full_Model//X_test_Morgan.pkl", "rb") as f:
    X_test = pickle.load(f)

with open("./path/to/intermediate/data/3.Constract_Full_Model/y_train_Morgan.pkl", "rb") as f:
    y_train = pickle.load(f)

with open("./path/to/intermediate/data/3.Constract_Full_Model/y_test_Morgan.pkl", "rb") as f:
    y_test = pickle.load(f)

# 3. Model construction (logistic regression / Random Forest)
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)

# 4. Prediction and evaluation
y_pred = rf.predict(X_test)

print("=== Classification results ===")
print(classification_report(y_test, y_pred))
print("=== Confusion matrix ===")
print(confusion_matrix(y_test, y_pred))

# === Classification results ===
# >>> print(classification_report(y_test, y_pred))
#               precision    recall  f1-score   support
# 
#            1       1.00      0.75      0.86         4
#            2       0.88      1.00      0.93         7
# 
#     accuracy                           0.91        11
#    macro avg       0.94      0.88      0.90        11
# weighted avg       0.92      0.91      0.91        11
# 
# >>> print("=== Confusion matrix ===")
# === Confusion matrix ===
# >>> print(confusion_matrix(y_test, y_pred))
# [[3 1]
#  [0 7]]

# Define output file name
with open("./path/to/intermediate/data/3.Constract_Full_Model/rf_model_Morgan.pkl", "wb") as f:
    pickle.dump(rf, f)

