##### ML FOR BIOSURVEILLANCE
# Perform Monte Carlo simulation-based classification of environments using metagenomic species-level features.

import pandas as pd

# Load input files
env_meta_data_file_path = "/content/Metadata_Aligned_to_species_corrected.csv"
env_data_file_path = "/content/Env_Pathogenic_species.csv"

# Load species abundance data
env_data = pd.read_csv(env_data_file_path)

# Load metadata with group and sample info
env_meta_data = pd.read_csv(env_meta_data_file_path)

# Preview data structure
print(env_data.head())
print(env_meta_data.head())

# -----------------------------------------------
# Function to Reshape Wide-Format to Long Format
# -----------------------------------------------

def combine_data(data, meta_data, columns_to_drop=[]):
    combine_df = []
    meta_data = meta_data.drop(columns=columns_to_drop)

    for i in range(len(meta_data)):
        row = meta_data.iloc[i]
        sample_id = row['Sample_ID']
        df = data[["taxon", sample_id]].rename(columns={sample_id: "relative_ab"})

        # Repeat metadata row to match species count
        repeated_row = pd.concat([row.to_frame().T]*len(df), ignore_index=True)
        df_id = pd.concat([df, repeated_row], axis=1)

        # Filter out species with 0 abundance
        df_id = df_id[df_id["relative_ab"] != 0.0]
        combine_df.append(df_id)

    return pd.concat(combine_df)

# Drop non-essential columns
columns_to_drop = [
    "Instrument", "Project_ID", "Sequencing_Center",
    "environmental_material", "Region", "Continent", "Country"
]

# Combine metadata and species matrix into long-form
Env_combine_df = combine_data(env_data, env_meta_data, columns_to_drop=columns_to_drop)
print(Env_combine_df.head())

# -----------------------------------------------
# Create Classification Matrix: Features and Labels
# -----------------------------------------------

from sklearn.preprocessing import LabelEncoder

# Pivot long-form data to wide format: samples × species
pivot_df = Env_combine_df.pivot_table(
    index=["Sample_ID", "Group"],
    columns="taxon",
    values="relative_ab",
    aggfunc="sum",
    fill_value=0
).reset_index()

# Extract X and y
X = pivot_df.drop(columns=["Sample_ID", "Group"])
y = pivot_df["Group"]

# Encode categorical features (if any remain)
for col in X.columns:
    if X[col].dtype == "object":
        X[col] = LabelEncoder().fit_transform(X[col])

print(X.shape)
print(y.value_counts())

# -----------------------------------------------
# Monte Carlo Classification with Random Forest
# -----------------------------------------------

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, cohen_kappa_score, confusion_matrix
from tqdm import tqdm
import numpy as np

# Encode labels
le = LabelEncoder()
y_encoded = le.fit_transform(y)
class_names = le.classes_

# Settings for simulation
n_iterations = 100
test_size = 0.2
metrics = []
conf_matrices = []

for i in tqdm(range(n_iterations), desc="Monte Carlo Simulation"):
    # Train/test split with stratification
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_encoded, stratify=y_encoded, test_size=test_size, random_state=i
    )

    # Train random forest
    model = RandomForestClassifier(
        class_weight='balanced', n_estimators=100,
        max_features='log2', min_samples_leaf=1,
        min_samples_split=2, max_depth=None, random_state=i
    )
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    # Store metrics
    acc = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, average="macro")
    kappa = cohen_kappa_score(y_test, y_pred)
    cm = confusion_matrix(y_test, y_pred, labels=range(len(class_names)))

    metrics.append([acc, f1, kappa])
    conf_matrices.append(cm)

# Save simulation results
metrics_df = pd.DataFrame(metrics, columns=["Accuracy", "F1_Macro", "Kappa"])
metrics_df.to_csv("monte_carlo_metrics.csv", index=False)  # Suppl. Table 7

# Average confusion matrix across simulations
conf_stack = np.stack(conf_matrices)
avg_cm = conf_stack.mean(axis=0)
avg_cm_norm = avg_cm / avg_cm.sum(axis=1, keepdims=True)

avg_cm_df = pd.DataFrame(avg_cm_norm, index=class_names, columns=class_names)
avg_cm_df.to_csv("average_confusion_matrix.csv")

# -----------------------------------------------
# Plot Confusion Matrix – Supplementary Figure 5
# -----------------------------------------------

import seaborn as sns
import matplotlib.pyplot as plt

# Load normalized confusion matrix
avg_cm_df = pd.read_csv("average_confusion_matrix.csv", index_col=0)

plt.figure(figsize=(8, 5))
sns.set_theme(style="whitegrid")

# Heatmap visualization
ax = sns.heatmap(
    avg_cm_df, annot=True, fmt=".2f", cmap="rocket_r",
    linewidths=0.5, linecolor="white", square=True,
    cbar_kws={"shrink": 0.85, "label": "Normalized Prediction Frequency"},
    annot_kws={"size": 5, "weight": "bold", "color": "black"}
)

# Customize labels and layout
plt.xticks(rotation=45, ha="right", fontsize=10, color="black")
plt.yticks(rotation=0, fontsize=10, color="black")
plt.xlabel("Predicted Environment", fontsize=10, color="black")
plt.ylabel("True Environment", fontsize=10, color="black")

# Save figure
plt.tight_layout()
plt.savefig("Confusion_Matrix(Figure 2).tiff", dpi=600, bbox_inches="tight")  # Suppl. Fig 5
plt.show()


########################
### Run Monte Carlo simulations to evaluate species-level feature importance 
########for classifying environments. Visualize stable predictors and ecological fingerprints.
###########################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from tqdm import tqdm
from sklearn.preprocessing import LabelEncoder

# ----------------------------------------
# STEP 1: Monte Carlo Simulation of Feature Importance
# ----------------------------------------

n_iterations = 100
all_importances = []

for i in tqdm(range(n_iterations), desc="Monte Carlo with Feature Importance"):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_encoded, stratify=y_encoded, test_size=0.2, random_state=i
    )

    model = RandomForestClassifier(
        n_estimators=100,
        max_features='log2',
        max_depth=None,
        min_samples_leaf=1,
        min_samples_split=2,
        class_weight='balanced',
        random_state=i
    )
    model.fit(X_train, y_train)

    # Store feature importances
    importances = model.feature_importances_
    all_importances.append(importances)

# ----------------------------------------
# STEP 2: Convert and Save Feature Importances
# ----------------------------------------

importance_array = np.array(all_importances)  # shape = [100 iterations × n_features]
species_names = X.columns

# Convert to DataFrame
importance_df = pd.DataFrame(importance_array, columns=species_names)

# Compute mean and standard deviation across iterations
summary_df = pd.DataFrame({
    "Species": species_names,
    "Mean_Importance": importance_df.mean(),
    "Std_Importance": importance_df.std()
}).sort_values(by="Mean_Importance", ascending=False)

summary_df.to_csv("monte_carlo_species_importance.csv", index=False)  # Supplementary Table 7(monte_carlo_species_importance)

# ----------------------------------------
# STEP 3: Bar Plot – Top 50 Most Predictive Species
# ----------------------------------------

summary_df = pd.read_csv("monte_carlo_species_importance.csv")
top_n = 50
top_df = summary_df.head(top_n)

sns.set_theme(style="white")
plt.figure(figsize=(9, 7))

bars = plt.barh(
    y=top_df["Species"],
    width=top_df["Mean_Importance"],
    xerr=top_df["Std_Importance"],
    color="#87CEEB", edgecolor="black", linewidth=0.7, capsize=2
)
plt.gca().invert_yaxis()
plt.xlabel("Mean Importance (± Std Dev)", fontsize=11, color='black')
plt.xticks(fontsize=10, color='black')
plt.yticks(fontsize=9, color='black')

# Beautify borders
for spine in plt.gca().spines.values():
    spine.set_visible(True)
    spine.set_linewidth(2)

plt.tight_layout()
plt.savefig("General_FeatureImportance.tiff", dpi=600, bbox_inches="tight", format="tiff")  # Suppl. Fig 6
plt.show()

# ----------------------------------------
# STEP 4: Boxplot of Feature Importance Distributions
# ----------------------------------------

# Convert to long-form
long_df = importance_df.melt(var_name="Species", value_name="Importance")
long_df = long_df[long_df["Species"].isin(summary_df.head(top_n)["Species"])]

plt.figure(figsize=(9, 7))
ax = sns.boxplot(
    data=long_df,
    y="Species",
    x="Importance",
    color="skyblue",
    linewidth=1.2,
    fliersize=2,
    orient="h"
)
plt.gca().invert_yaxis()
plt.xlabel("Feature Importance", fontsize=10, color="black")
plt.ylabel("", fontsize=10)
plt.xticks(fontsize=10, color="black")
plt.yticks(fontsize=9, color="black")
plt.grid(False)

# Thicker border
for spine in plt.gca().spines.values():
    spine.set_visible(True)
    spine.set_linewidth(2)

plt.tight_layout()
plt.savefig("Feature_Importance_Distribution(Suppl 4).tiff", dpi=600, format="tiff", bbox_inches="tight")
plt.show()

# ----------------------------------------
# STEP 5: Ecological Fingerprint Heatmap – Top 100 Predictive Species
# ----------------------------------------

# Drop non-numeric metadata columns
features_only = pivot_df.drop(columns=["Sample_ID", "Group"])
features_only["Group"] = pivot_df["Group"]

# Compute mean relative abundance of species by environment
species_means = features_only.groupby("Group").mean()

# Normalize by z-score per species
species_fingerprint = species_means.T.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# Filter for top 100 species
top_species = summary_df.head(100)["Species"]
species_fingerprint_top = species_fingerprint.loc[top_species]

# Plot heatmap
plt.figure(figsize=(10, 12))
sns.set_theme(style="white")
ax = sns.heatmap(
    species_fingerprint_top,
    cmap="vlag", center=0,
    linewidths=0.4, linecolor="gray",
    cbar_kws={"label": "Z-score (Normalized Abundance)"},
    xticklabels=True, yticklabels=True
)
plt.xlabel("", fontsize=13)
plt.ylabel("", fontsize=13)
plt.xticks(fontsize=10, rotation=30, color="black")
plt.yticks(fontsize=8, color="black")

# Thicker frame
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(2)

plt.tight_layout()
plt.savefig("Ecological_Fingerprint_Heatmap_100.tiff", dpi=600, format="tiff", bbox_inches="tight")  # Figure 4
plt.show()

########################################
#### Evaluate microbial resilience, contamination risk, and prediction uncertainty using synthetic microbial mixtures
#########################################
# Import standard libraries
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import math
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import LabelEncoder
from scipy.stats import entropy

# ----------------------------
# STEP 1: Synthetic Microbial Mixing Simulation
# ----------------------------

# Assumes pivot_df, X, y_encoded, le, model are already defined
env_names = pivot_df["Group"].unique()
species_cols = [col for col in pivot_df.columns if col not in ["Sample_ID", "Group"]]
grouped = pivot_df.groupby("Group")

blend_ratios = np.linspace(0, 1, 11)  # 0% to 100% mixing in 10% steps
results = []

# Loop through all environment pairs and simulate blends
for env_a, env_b in itertools.combinations(env_names, 2):
    print(f"Simulating: {env_a} ↔ {env_b}")
    samples_a = grouped.get_group(env_a)[species_cols].values
    samples_b = grouped.get_group(env_b)[species_cols].values
    n_samples = min(len(samples_a), len(samples_b))

    for ratio in blend_ratios:
        for i in range(n_samples):
            blend = ratio * samples_a[i] + (1 - ratio) * samples_b[i]
            proba = model.predict_proba([blend])[0]
            predicted_class = np.argmax(proba)
            predicted_label = le.inverse_transform([predicted_class])[0]
            confidence = np.max(proba)
            ent = entropy(proba, base=2)

            results.append({
                "Env_A": env_a,
                "Env_B": env_b,
                "Blend_A": ratio,
                "Blend_B": 1 - ratio,
                "Predicted": predicted_label,
                "Confidence": confidence,
                "Entropy": ent
            })

# Save blending results
blend_df = pd.DataFrame(results)
blend_df.to_csv("synthetic_blend_entropy.csv", index=False)  # Supplementary Table 9 (synthetic_blend_entropy)


# ----------------------------
# STEP 2: Plot Confidence Curves – Figure 5
# ----------------------------

blend_df = pd.read_csv("synthetic_blend_entropy.csv")

# Setup subplots
env_pairs = blend_df[["Env_A", "Env_B"]].drop_duplicates()
n_plots = len(env_pairs)
cols = 3
rows = math.ceil(n_plots / cols)
fig, axes = plt.subplots(rows, cols, figsize=(cols * 6.5, rows * 5), squeeze=False)
palette = sns.color_palette("Set2")

# Plot each environment pair's confidence curves
for idx, ((env_a, env_b), subset) in enumerate(blend_df.groupby(["Env_A", "Env_B"])):
    row, col = idx // cols, idx % cols
    ax = axes[row][col]

    sns.lineplot(data=subset, x="Blend_A", y="Confidence", hue="Predicted",
                 palette=palette, marker="o", ax=ax, linewidth=2.2, alpha=0.9)

    ax.axvline(0.5, linestyle="--", color="gray", lw=1, alpha=0.7)
    ax.set_title(f"{env_a} ↔ {env_b}", fontsize=13, weight="bold")
    ax.set_xlabel(f"% {env_a}", fontsize=11, color="black")
    ax.set_ylabel("Prediction Confidence", fontsize=11, color="black")
    ax.set_ylim(0, 1.05)
    ax.legend(loc='lower left', fontsize=9)

# Remove unused subplots
for j in range(idx + 1, rows * cols):
    fig.delaxes(axes[j // cols][j % cols])

plt.tight_layout(rect=[0, 0, 1, 0.98])
plt.savefig("Synthetic_Mixing_Curves.tiff", dpi=600, format="tiff", bbox_inches="tight")  # Figure 5
plt.show()



# ----------------------------
# STEP 3: Extract Switch Points 
# ----------------------------

# Identify first switch point per blending simulation
switch_points = []

for (env_a, env_b), group in blend_df.groupby(["Env_A", "Env_B"]):
    group_sorted = group.sort_values("Blend_A")
    predictions = group_sorted["Predicted"].values
    blend_percents = group_sorted["Blend_A"].values

    for i in range(1, len(predictions)):
        if predictions[i] != predictions[i - 1]:
            switch_points.append({
                "Env_A": env_a,
                "Env_B": env_b,
                "Switch_From": predictions[i - 1],
                "Switch_To": predictions[i],
                "Blend_Percent_A": blend_percents[i],
                "Blend_Percent_B": 1 - blend_percents[i]
            })
            break

switch_df = pd.DataFrame(switch_points)
switch_df.to_csv("synthetic_transition_points.csv", index=False)  # Supplementary Table 8 (synthetic_transition_points)


# ----------------------------
# STEP 4: Plot Blend Switch Heatmap 
# ----------------------------

switch_df = pd.read_csv("synthetic_transition_points.csv")

# Build symmetric matrix of switch blend ratios
envs = sorted(set(switch_df["Env_A"]) | set(switch_df["Env_B"]))
matrix = pd.DataFrame(index=envs, columns=envs, dtype=float)

for _, row in switch_df.iterrows():
    a, b = row["Env_A"], row["Env_B"]
    matrix.loc[a, b] = row["Blend_Percent_A"]
    matrix.loc[b, a] = row["Blend_Percent_B"]

plt.figure(figsize=(7, 5))
sns.set_theme(style="white")
ax = sns.heatmap(
    matrix.astype(float),
    annot=True,
    cmap="vlag",
    fmt=".2f",
    linewidths=0.5,
    linecolor="white",
    cbar_kws={"label": "Blend Ratio", "shrink": 0.8}
)
plt.xlabel("To Environment", fontsize=12, weight="bold", color="black")
plt.ylabel("From Environment", fontsize=12, weight="bold", color="black")
plt.tight_layout()
plt.savefig("Blend_Switch_Point_Heatmap(Figure 6).tiff", dpi=600, format="tiff", bbox_inches="tight")  ### Figure 6A
plt.show()



# ----------------------------
# STEP 5: Dominance Matrix – Figure 6B
# ----------------------------

dominance_matrix = pd.DataFrame(0, index=envs, columns=envs)

for _, row in switch_df.iterrows():
    winner, loser = row["Switch_To"], row["Switch_From"]
    dominance_matrix.loc[winner, loser] += 1

dominance_matrix.to_csv("environment_dominance_matrix.csv")

plt.figure(figsize=(8, 6))
sns.heatmap(
    dominance_matrix,
    annot=True,
    fmt="d",
    cmap="YlGnBu",
    linewidths=0.5,
    cbar_kws={"label": "Number of Wins", "shrink": 0.85}
)
plt.xlabel("Loser", fontsize=12, weight="bold")
plt.ylabel("Winner", fontsize=12, weight="bold")
plt.tight_layout()
plt.savefig("Environment_Dominance_Heatmap(Figure 7).tiff", dpi=600, format="tiff", bbox_inches="tight")
plt.show()

# ----------------------------
# STEP 6: Entropy Curves – Figure 8
# ----------------------------

blend_df = pd.read_csv("synthetic_blend_entropy.csv")

env_pairs = blend_df[["Env_A", "Env_B"]].drop_duplicates()
n_plots = len(env_pairs)
cols, rows = 3, math.ceil(n_plots / 3)
fig, axes = plt.subplots(rows, cols, figsize=(cols * 6.5, rows * 4.5), squeeze=False)

max_entropy = np.log2(blend_df["Predicted"].nunique())

for idx, ((env_a, env_b), subset) in enumerate(blend_df.groupby(["Env_A", "Env_B"])):
    row, col = idx // cols, idx % cols
    ax = axes[row][col]

    summary = subset.groupby("Blend_A")["Entropy"].agg(["mean", "std"]).reset_index()

    sns.lineplot(data=summary, x="Blend_A", y="mean", ax=ax, linewidth=2, marker="o")
    ax.fill_between(summary["Blend_A"], summary["mean"] - summary["std"], summary["mean"] + summary["std"], alpha=0.25, color="steelblue")
    ax.axhline(max_entropy, linestyle="--", color="red", linewidth=1.5)
    
    # Switch point
    sw = switch_df[(switch_df["Env_A"] == env_a) & (switch_df["Env_B"] == env_b)]
    if not sw.empty:
        ax.axvline(sw.iloc[0]["Blend_Percent_A"], linestyle="--", color="black", linewidth=1.2)

    ax.set_title(f"{env_a} ↔ {env_b}", fontsize=12, weight="bold")
    ax.set_xlabel(f"% {env_a}", fontsize=12)
    ax.set_ylabel("Entropy", fontsize=12)
    ax.set_ylim(0, max_entropy + 0.2)

for j in range(idx + 1, rows * cols):
    fig.delaxes(axes[j // cols][j % cols])

plt.tight_layout(rect=[0, 0, 1, 0.98])
plt.savefig("Synthetic_Blend_Entropy_Curves(Figure 8).tiff", dpi=600, bbox_inches="tight", format="tiff")
plt.show()

####################
### Data-driven contamination risk index by combining multiple metrics
####################

# ============================
# STEP 1: Load Precomputed Inputs
# ============================

import pandas as pd

# Load switch matrix for environment transition dominance
# (already created from synthetic blending switch point analysis)
switch_df = pd.read_csv("synthetic_transition_points.csv")
blend_df = pd.read_csv("synthetic_blend_entropy.csv")

# Count how many times each environment was overridden (i.e., lost prediction)
loser_counts = switch_df["Switch_From"].value_counts()
loser_counts.name = "Times_Overridden"

# Compute mean entropy for each environment across all blends it participated in
env_entropy = {}
for env in blend_df["Env_A"].unique():
    mask = (blend_df["Env_A"] == env) | (blend_df["Env_B"] == env)
    mean_entropy = blend_df[mask].groupby("Blend_A")["Entropy"].mean().mean()
    env_entropy[env] = mean_entropy

entropy_df = pd.DataFrame.from_dict(env_entropy, orient="index", columns=["Mean_Entropy"])

# Combine both into a single DataFrame
risk_df = pd.concat([loser_counts, entropy_df], axis=1).fillna(0)


# ============================
# STEP 2: Normalize + Basic Contamination Risk Score
# ============================

# Normalize each metric (min-max scaling)
risk_df["Norm_Overridden"] = (risk_df["Times_Overridden"] - risk_df["Times_Overridden"].min()) / (risk_df["Times_Overridden"].max() - risk_df["Times_Overridden"].min())
risk_df["Norm_Entropy"] = (risk_df["Mean_Entropy"] - risk_df["Mean_Entropy"].min()) / (risk_df["Mean_Entropy"].max() - risk_df["Mean_Entropy"].min())

# Weighted risk score (tunable weights)
risk_df["Contamination_Risk_Score"] = 0.6 * risk_df["Norm_Overridden"] + 0.4 * risk_df["Norm_Entropy"]

# Save sorted table
risk_df = risk_df.sort_values("Contamination_Risk_Score", ascending=False)
risk_df.to_csv("environment_contamination_risk.csv")    


# ============================
# STEP 3: Barplot – Risk Score by Environment 
# ============================

import matplotlib.pyplot as plt
import seaborn as sns

# Load for plotting
risk_df = pd.read_csv("environment_contamination_risk.csv", index_col=0)

plt.figure(figsize=(5, 2))
sns.set_theme(style="white")

ax = sns.barplot(
    data=risk_df.reset_index(),
    x="Contamination_Risk_Score",
    y="index",
    palette="Reds_r",
    edgecolor="black"
)

# Format labels
plt.xlabel("Risk Score", fontsize=8, color="black")
plt.ylabel("", fontsize=10, color="black")
plt.xticks(fontsize=8, color="black")
plt.yticks(fontsize=8, color="black")
plt.xlim(0, 1)

# Make borders thicker
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(2)

plt.grid(False)
plt.tight_layout()
plt.savefig("Contamination_Risk_Barplot.tiff", dpi=600, format="tiff", bbox_inches="tight")   ### Figure 6D
plt.show()



# ============================
# STEP 4: Add Richness and Variance to Risk Model
# ============================

richness_dict = {}
variance_dict = {}

for env in pivot_df["Group"].unique():
    samples = pivot_df[pivot_df["Group"] == env][species_cols]
    richness = (samples > 0).sum(axis=1).mean()  # Avg # of species per sample
    variance = samples.var(axis=1).mean()        # Mean variance across samples
    richness_dict[env] = richness
    variance_dict[env] = variance

# Create new columns and merge into risk_df
richness_df = pd.DataFrame.from_dict(richness_dict, orient="index", columns=["Richness"])
variance_df = pd.DataFrame.from_dict(variance_dict, orient="index", columns=["Variance"])
risk_df = pd.concat([risk_df, richness_df, variance_df], axis=1).fillna(0)

# Normalize all columns
for col in ["Times_Overridden", "Mean_Entropy", "Richness", "Variance"]:
    risk_df[f"Norm_{col}"] = (risk_df[col] - risk_df[col].min()) / (risk_df[col].max() - risk_df[col].min())

# Extended scoring models
risk_df["Contamination_Risk_Score_50_50"] = 0.5 * risk_df["Norm_Times_Overridden"] + 0.5 * risk_df["Norm_Mean_Entropy"]

risk_df["Contamination_Risk_Score_Expanded"] = (
    0.4 * risk_df["Norm_Times_Overridden"] +
    0.3 * risk_df["Norm_Mean_Entropy"] +
    0.15 * risk_df["Norm_Richness"] +
    0.15 * risk_df["Norm_Variance"]
)

# Save final extended table
risk_df = risk_df.sort_values("Contamination_Risk_Score_Expanded", ascending=False)
risk_df.to_csv("extended_environment_risk_profile.csv")  # Supplementary Table 9 (Environemnt Risk Profile)


# ============================
# STEP 5: Radar Plot of Risk Dimensions 
# ============================

import numpy as np
import matplotlib.pyplot as plt

risk_df = pd.read_csv("extended_environment_risk_profile.csv", index_col=0)

# Select the 4 dimensions
categories = ["Overridden", "Entropy", "Richness", "Variance"]
radar_data = risk_df[["Norm_Times_Overridden", "Norm_Mean_Entropy", "Norm_Richness", "Norm_Variance"]]
radar_data.columns = categories
radar_data = radar_data.reset_index()

# Setup radar chart
num_vars = len(categories)
angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
angles += angles[:1]

colors = ['#E76F51', '#2A9D8F', '#264653', '#F4A261']

plt.figure(figsize=(10, 8))
ax = plt.subplot(111, polar=True)
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

# Plot each environment
for i, row in radar_data.iterrows():
    values = row[categories].tolist()
    values += values[:1]
    label = row["index"]

    ax.plot(angles, values, color=colors[i], linewidth=2)
    ax.fill(angles, values, color=colors[i], alpha=0.25)

    # Add text label inline
    angle = angles[i % len(categories)]
    radius = values[i % len(categories)]
    ax.text(angle, radius + 0.02, label, fontsize=11, weight='bold', color=colors[i], ha='center')

# Styling
ax.set_ylim(0, 1.0)
ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories, fontsize=12, weight='bold')
ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=10)
ax.yaxis.grid(True, linestyle='dotted')
ax.xaxis.grid(True, linestyle='--')
ax.spines["polar"].set_visible(False)

# Outer border
plt.gca().patch.set_edgecolor('black')
plt.gca().patch.set_linewidth(2)

plt.tight_layout()
plt.savefig("Ecological_Risk_Radar_Plot.tiff", dpi=600, format="tiff", bbox_inches="tight")  # Figure 6C
plt.show()







































