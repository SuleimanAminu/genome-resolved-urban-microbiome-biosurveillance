import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ttest_ind, entropy

# -----------------------------------------------
# STEP 1: Load input files
# -----------------------------------------------

# Load wide-format pathogenic species matrix and aligned metadata
species_matrix_path = "your_path/Env_Pathogenic_species.csv"
metadata_path = "your_path/Metadata_Aligned_to_species_corrected.csv"

species_df = pd.read_csv(species_matrix_path)
metadata_df = pd.read_csv(metadata_path)

# Rename first column to "taxon" (species name)
species_df.rename(columns={species_df.columns[0]: "taxon"}, inplace=True)

# Drop rows without species name
species_df_clean = species_df.dropna(subset=["taxon"])

# -----------------------------------------------
# STEP 2: Convert to long format and merge metadata
# -----------------------------------------------

# Melt species abundance matrix to long format
species_long_fixed = species_df_clean.melt(id_vars=["taxon"], var_name="Sample_ID", value_name="Abundance")

# Merge with metadata to get environmental Group info
species_long_fixed = species_long_fixed.merge(metadata_df, on="Sample_ID")

# -----------------------------------------------
# STEP 3: Calculate prevalence per species per group
# -----------------------------------------------

prevalence_df_corrected = species_long_fixed.groupby(["taxon", "Group"]).apply(
    lambda x: (x["Abundance"] > 0).sum() / len(x) * 100
).reset_index(name="Prevalence (%)")

# -----------------------------------------------
# STEP 4: Define key pathogens (exact + grouped)
# -----------------------------------------------

exact_matches = [
    "Escherichia coli", "Acinetobacter baumannii", "Mycobacterium tuberculosis", "Salmonella Typhii",
    "Enterococcus faecium", "Pseudomonas aeruginosa", "Neisseria gonorrhoeae",
    "Staphylococcus aureus", "Streptococcus pneumoniae", "Haemophilus influenza",
    "Helicobacter pylori", "Clostridioides difficile", "Klebsiella oxytoca", 
    "Enterococcus faecalis", "Enterococcus avium", "Candida albicans", 
    "Candida parapsilosis", "Candida glabrata", "Candida dubliniensis", 
    "Proteus mirabilis"
]

group_matches = {
    "Shigella species": "Shigella",
    "Enterobacter species": "Enterobacter",
    "Citrobacter species": "Citrobacter",
    "Proteus species": "Proteus",
    "Serratia species": "Serratia",
    "Streptococcus species": "Streptococcus",
    "Morganella species": "Morganella",
    "Providencia species": "Providencia",
    "Campylobacter species": "Campylobacter",
}

# -----------------------------------------------
# STEP 5: Filter and group prevalence data
# -----------------------------------------------

# Filter exact matches
filtered_exact = prevalence_df_corrected[prevalence_df_corrected["taxon"].isin(exact_matches)]

# Filter species that match grouped patterns
filtered_groups = prevalence_df_corrected[
    prevalence_df_corrected["taxon"].apply(lambda x: any(x.startswith(g) for g in group_matches.values()))
].copy()

# Replace species names with their grouped label
for group_name, prefix in group_matches.items():
    filtered_groups.loc[filtered_groups["taxon"].str.startswith(prefix), "taxon"] = group_name

# Take max prevalence per group/species
filtered_groups = filtered_groups.groupby(["taxon", "Group"], as_index=False).agg({"Prevalence (%)": "max"})

# Combine exact and grouped
final_filtered_df = pd.concat([filtered_exact, filtered_groups])

# Include manually added species if missing
additional_species = ["Klebsiella pneumoniae", "Stenotrophomonas maltophilia"]
filtered_additional = prevalence_df_corrected[prevalence_df_corrected["taxon"].isin(additional_species)]
final_filtered_df_updated = pd.concat([final_filtered_df, filtered_additional])

# Select top 30 most prevalent species overall
top_final_species = final_filtered_df_updated.groupby("taxon")["Prevalence (%)"].max().nlargest(30).index
filtered_final_top = final_filtered_df_updated[final_filtered_df_updated["taxon"].isin(top_final_species)]

# Save results
filtered_final_top.to_csv("Group_Pathogenic_Species_Prevalence.csv", index=False)  # Supplementary Table 5

# -----------------------------------------------
# STEP 6: Plot stacked horizontal bar chart – Fig 3E
# -----------------------------------------------

# Reformat to wide matrix for plotting
stacked_data = filtered_final_top.pivot(index="taxon", columns="Group", values="Prevalence (%)").fillna(0)
species_names = stacked_data.index[::-1]
y_pos = np.arange(len(species_names))

colors = ["#D73027", "#FC8D59", "#FEE090", "#91CF60", "#1A9850", "#3288BD", "#BEBADA", "#FDB462"]

fig, ax = plt.subplots(figsize=(5, 3))
bottom_values = np.zeros(len(species_names))

for i, group in enumerate(stacked_data.columns):
    ax.barh(
        y_pos + 0.15, stacked_data.loc[species_names, group], left=bottom_values,
        label=group, color=colors[i % len(colors)], alpha=0.85,
        edgecolor="black", linewidth=0.8, height=0.8
    )
    bottom_values += stacked_data.loc[species_names, group]

# Format bar plot
ax.set_yticks(y_pos + 0.15)
ax.set_yticklabels(species_names, fontsize=7, color="black")
ax.set_xlabel("Prevalence (%)", fontsize=8, color="black")
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
plt.setp(ax.get_xticklabels(), fontsize=8, color="black")
ax.legend(loc="upper right", fontsize=5, frameon=False, bbox_to_anchor=(1, 0.7))
ax.invert_yaxis()

# Clean aesthetics
for spine in ax.spines.values():
    spine.set_color("black")
    spine.set_linewidth(2)

plt.grid(False)
plt.savefig("WHO_Group_Prevalence.png", dpi=600, bbox_inches="tight")  # Figure 3E
plt.show()

# -----------------------------------------------
# STEP 7: Violin plot of pathogen prevalence – Fig 3D
# -----------------------------------------------

fig, ax = plt.subplots(figsize=(2, 3))
colors = sns.color_palette("tab20", 17)
sns.violinplot(data=filtered_final_top, x="Group", y="Prevalence (%)", palette=colors, ax=ax)

plt.xticks(rotation=45, ha='right', fontsize=10, color="black")
plt.xlabel("", fontsize=2)
plt.ylabel("Prevalence (%)", fontsize=8, color="black")

# Border styling
for spine in ax.spines.values():
    spine.set_color("black")
    spine.set_linewidth(2)

plt.grid(False)
plt.savefig("WHO_Group_General_Prevalence.png", dpi=600, bbox_inches="tight")  # Figure 3D
plt.show()

# -----------------------------------------------
# STEP 8: Compute Shannon diversity – For comparison
# -----------------------------------------------

# Create sample × species matrix
abundance_matrix = species_long_fixed.pivot_table(index='Sample_ID', columns='taxon', values='Abundance', fill_value=0)

# Calculate Shannon diversity
abundance_matrix['Shannon_Diversity'] = abundance_matrix.apply(lambda row: entropy(row + 1e-9, base=np.e), axis=1)

# Save per-sample Shannon
sample_diversity = abundance_matrix[['Shannon_Diversity']].reset_index()

# -----------------------------------------------
# STEP 9: Compare diversity based on presence of key pathogens – Supplementary Figure 4
# -----------------------------------------------

target_pathogens = [
    "Escherichia coli", "Enterococcus faecium", "Staphylococcus aureus",
    "Klebsiella pneumoniae", "Stenotrophomonas maltophilia",
    "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Enterococcus faecalis"
]

diversity_results = {}
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(16, 10))
axes = axes.flatten()

for i, pathogen in enumerate(target_pathogens):
    species_long_fixed["Pathogen_Presence"] = (species_long_fixed["taxon"] == pathogen) & (species_long_fixed["Abundance"] > 0)
    pathogen_presence_df = species_long_fixed.groupby("Sample_ID")["Pathogen_Presence"].any().reset_index(name="Pathogen_Present")
    sample_diversity_with_pathogen = sample_diversity.merge(pathogen_presence_df, on="Sample_ID")

    shannon_present = sample_diversity_with_pathogen[sample_diversity_with_pathogen["Pathogen_Present"]]["Shannon_Diversity"]
    shannon_absent = sample_diversity_with_pathogen[~sample_diversity_with_pathogen["Pathogen_Present"]]["Shannon_Diversity"]

    t_test_result = ttest_ind(shannon_present, shannon_absent, equal_var=False, nan_policy='omit')

    diversity_results[pathogen] = {
        "Mean Diversity (Present)": shannon_present.mean(),
        "Mean Diversity (Absent)": shannon_absent.mean(),
        "P-Value": t_test_result.pvalue
    }

    sns.boxplot(x="Pathogen_Present", y="Shannon_Diversity", data=sample_diversity_with_pathogen, ax=axes[i])
    axes[i].set_xlabel("Pathogen Presence")
    axes[i].set_ylabel("Shannon Diversity")
    axes[i].set_title(f"{pathogen} (p={t_test_result.pvalue:.4f})")

plt.tight_layout()
plt.savefig("Alpha_Diversity_Pathogen_Comparison.png", dpi=300)  # Supplementary Figure 4
plt.show()

diversity_results_df = pd.DataFrame.from_dict(diversity_results, orient="index")
diversity_results_df.to_csv("Diversity_Result_Pathogen_influence.csv", index=False)  # Supplementary Table 5
print(diversity_results_df)

