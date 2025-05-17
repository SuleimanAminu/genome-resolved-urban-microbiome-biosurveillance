# -----------------------------------------------
# Script: 02_model_comparison.py
# Authors: Suleiman & AbdulAziz
# Description: Compare multiple classifiers on environment classification task
# -----------------------------------------------

import pandas as pd
import numpy as np
import os
from sklearn.model_selection import cross_validate, StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import make_scorer, cohen_kappa_score
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC, LinearSVC
from sklearn.neighbors import KNeighborsClassifier

# === Load Data ===
base_path = "/srv/lustre01/project/mmrd-cp3fk69sfrq/user"
abundance_file = os.path.join(base_path, "Env_Pathogenic_species.csv")
metadata_file = os.path.join(base_path, "Metadata_Aligned_to_species_corrected.csv")

print("Loading data...")
env_data = pd.read_csv(abundance_file)
env_meta = pd.read_csv(metadata_file)

# === Combine metadata and abundance ===
def combine_data(data, meta_data, columns_to_drop=[]):
    combined_list = []
    meta_data = meta_data.drop(columns=columns_to_drop)
    for i in range(len(meta_data)):
        row = meta_data.iloc[i]
        sample_id = row["Sample_ID"]
        if sample_id not in data.columns:
            continue
        df = data[["Species", sample_id]].rename(columns={sample_id: "relative_ab"})
        repeated_row = pd.concat([row.to_frame().T] * len(df), ignore_index=True)
        df_id = pd.concat([df, repeated_row], axis=1)
        df_id = df_id[df_id["relative_ab"] != 0.0]
        combined_list.append(df_id)
    return pd.concat(combined_list, ignore_index=True).drop(columns=["Sample_ID"])

drop_cols = ["Instrument", "Project_ID", "Sequencing_Center",
             "environmental_material", "Region", "Continent", "Country"]

print(" Combining data...")
combined_df = combine_data(env_data, env_meta, columns_to_drop=drop_cols)
combined_df.to_csv(os.path.join(base_path, "combined_df.csv"), index=False)

# === Prepare features and labels ===
target = "Group"
X = combined_df.drop(columns=[target])
y = combined_df[target]

# Encode categorical features
for col in X.columns:
    if X[col].dtype == "object":
        X[col] = X[col].astype("category").cat.codes

# Encode target labels
label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y)

# === Define metrics ===
def kappa_scorer(y_true, y_pred):
    return cohen_kappa_score(y_true, y_pred)

scoring = {
    "accuracy": "accuracy",
    "f1_macro": "f1_macro",
    "kappa": make_scorer(kappa_scorer)
}

# === Stratified K-Fold ===
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

# === Define models ===
models = {
    "Regularized_RF": RandomForestClassifier(
        n_estimators=100, max_depth=None, max_features="log2",
        min_samples_leaf=1, min_samples_split=2,
        random_state=42, class_weight="balanced"
    ),
    "Standard_RF": RandomForestClassifier(random_state=42, class_weight="balanced"),
    "Logistic_Regression": LogisticRegression(max_iter=1000),
    "Linear_SVM": LinearSVC(max_iter=1000, random_state=42),
    "KNN": KNeighborsClassifier(n_neighbors=5),
    "Gradient_Boosting": GradientBoostingClassifier()
}

# === Cross-validation loop ===
results = {}
all_scores = []

for name, model in models.items():
    print(f" Evaluating: {name}")
    scores = cross_validate(model, X, y_encoded, cv=cv, scoring=scoring, n_jobs=-1)
    
    results[name] = {
        "Accuracy": np.mean(scores["test_accuracy"]),
        "F1 Macro": np.mean(scores["test_f1_macro"]),
        "Kappa": np.mean(scores["test_kappa"])
    }

    for i in range(len(scores["test_accuracy"])):
        all_scores.append({
            "Model": name,
            "Fold": i + 1,
            "Accuracy": scores["test_accuracy"][i],
            "F1 Macro": scores["test_f1_macro"][i],
            "Kappa": scores["test_kappa"][i]
        })

# === Save results ===
results_df = pd.DataFrame(results).T.sort_values(by="Kappa", ascending=False)
results_df.to_csv(os.path.join(base_path, "Env_model_comparison_results.csv"))

all_scores_df = pd.DataFrame(all_scores)
all_scores_df.to_csv(os.path.join(base_path, "Env_model_foldwise_scores.csv"), index=False)

# === Display summary ===
print("\n Model Comparison Summary (sorted by Kappa):")
print(results_df)

