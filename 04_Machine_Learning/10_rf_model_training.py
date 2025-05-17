# -----------------------------------------------
# Script: 01_rf_model_training.py
# Authors: Suleiman & AbdulAziz
# Description: Random Forest classifier training with GridSearchCV on environmental microbiome data
# -----------------------------------------------

import pandas as pd
import numpy as np
import joblib
import os
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.metrics import classification_report, confusion_matrix

# === Load Data ===
base_path = "/srv/lustre01/project/mmrd-cp3fk69sfrq/user"
meta_path = os.path.join(base_path, "Metadata_Aligned_to_species_corrected.csv")
abundance_path = os.path.join(base_path, "Env_Pathogenic_species.csv")

print(" Loading metadata and abundance data...")
metadata = pd.read_csv(meta_path)
abundance = pd.read_csv(abundance_path)

# === Combine Metadata and Species Data ===
def combine_data(data, meta_data, columns_to_drop=[]):
    combined = []
    meta_data = meta_data.drop(columns=columns_to_drop)
    for i in range(len(meta_data)):
        row = meta_data.iloc[i]
        sample_id = row['Sample_ID']
        if sample_id not in data.columns:
            continue
        df = data[["Species", sample_id]].rename(columns={sample_id: 'relative_ab'})
        repeated_row = pd.concat([row.to_frame().T]*len(df), ignore_index=True)
        df_combined = pd.concat([df, repeated_row], axis=1)
        df_combined = df_combined[df_combined['relative_ab'] != 0.0]
        combined.append(df_combined)
    return pd.concat(combined).drop(columns=['Sample_ID'])

columns_to_drop = ["Instrument", "Project_ID", "Sequencing_Center",
                   "environmental_material", "Region", "Continent", "Country"]

print("Combining datasets...")
df = combine_data(abundance, metadata, columns_to_drop=columns_to_drop)

# === Prepare Features and Labels ===
target = "Group"
X = df.drop(columns=[target])
y = df[target]

# Label encode object features (except Species, which we one-hot encode)
label_encoder = LabelEncoder()
for col in X.columns:
    if col != 'Species' and X[col].dtype == 'object':
        X[col] = label_encoder.fit_transform(X[col])

# One-hot encode Species
if 'Species' in X.columns:
    X = pd.get_dummies(X, columns=['Species'])

# === Train-Test Split ===
X_train, X_test, y_train, y_test = train_test_split(
    X, y, stratify=y, test_size=0.2, random_state=42
)

# === Define Parameter Grid ===
param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [None, 5, 10, 15, 20],
    'min_samples_split': [2, 5],
    'min_samples_leaf': [1, 2],
    'max_features': ['sqrt', 'log2']
}

scoring = {
    'accuracy': 'accuracy',
    'f1_macro': 'f1_macro'
}

print("Starting Random Forest training with GridSearchCV...")
rf = RandomForestClassifier(random_state=42, class_weight="balanced")
grid_search = GridSearchCV(
    estimator=rf,
    param_grid=param_grid,
    refit='f1_macro',
    scoring=scoring,
    cv=10,
    n_jobs=-1,
    verbose=2
)

grid_search.fit(X_train, y_train)
best_rf = grid_search.best_estimator_

print("Best Parameters:")
print(grid_search.best_params_)

# === Model Evaluation ===
y_pred = best_rf.predict(X_test)
print(" Classification Report:")
print(classification_report(y_test, y_pred))

conf_matrix = pd.DataFrame(confusion_matrix(y_test, y_pred),
                           index=np.unique(y_test),
                           columns=np.unique(y_test))

print("Confusion Matrix:")
print(conf_matrix)

# === Save the Trained Model ===
model_save_path = os.path.join(base_path, "best_random_forest_model.pkl")
joblib.dump(best_rf, model_save_path)
