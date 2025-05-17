import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#######################################


# Load datasets (update paths if needed)
file_paths = {
    "Public Transport": "your_path/Public_transport_selected_bin_cleaned.csv",
    "Hospital Environment": "your_path/Hosp_env_selected_bin_cleaned.csv",
    "Hospital Sewage": "your_path/Hosp_sewage_selected_bin_cleaned.csv",
    "PRJNA369713": "your_path/PRJNA369713_selected_bin_cleaned.csv",
}

# Read and combine all datasets
dfs = {name: pd.read_csv(path) for name, path in file_paths.items()}
combined_df = pd.concat(dfs.values(), ignore_index=True)

# Display dataset info
print("\nDataset Overview:\n")
print(combined_df.info())

# Display first few rows
print("\nFirst 5 Rows:\n", combined_df.head())

# Summary statistics
print("\nSummary Statistics:\n", combined_df.describe())


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Ensure numeric conversion (handling errors)
numeric_columns = ["Completeness", "Contamination", "Strain_heterogeneity"]
combined_df[numeric_columns] = combined_df[numeric_columns].apply(pd.to_numeric, errors="coerce")

# Drop any rows with NaN values in Completeness
cleaned_df = combined_df.dropna(subset=["Completeness"])

# Set plot style
sns.set_style("whitegrid")

# Create a figure with a smaller size and high resolution
fig, ax = plt.subplots(figsize=(5, 4))

# KDE Plot for Completeness Distribution by Category
sns.kdeplot(data=cleaned_df, x="Completeness", hue="Category", fill=True, alpha=0.5, ax=ax)

# Formatting the plot
#ax.set_title("Distribution of Completeness by Category", fontsize=12, fontweight='bold')
ax.set_xlabel("Completeness (%)", fontsize=10, color='black')
ax.set_ylabel("Density", fontsize=10)
ax.set_ylabel("Density", fontsize=10, color='black')
#ax.legend(title="Category", fontsize=8, title_fontsize=9)
#sns.despine()

# Remove grid lines
ax.grid(False)

# Make the borders black
for spine in ax.spines.values():
    spine.set_edgecolor("black")
    spine.set_linewidth(1.5)

# Save the figure at 600 DPI
save_path = "your_path/Env_MAGS_completeness_distribution.png"     ##### Figure 2D
plt.savefig(save_path, dpi=600, bbox_inches="tight")

# Show the plot
plt.show()

#####################################################
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Ensure numeric conversion (handling errors)
numeric_columns = ["Completeness", "Contamination", "Strain_heterogeneity"]
combined_df[numeric_columns] = combined_df[numeric_columns].apply(pd.to_numeric, errors="coerce")

# Drop any rows with NaN values in Completeness
cleaned_df = combined_df.dropna(subset=["Completeness"])

# Set plot style
sns.set_style("whitegrid")

# Create a figure with a smaller size and high resolution
fig, ax = plt.subplots(figsize=(5, 4))



# KDE Plot for Completeness Distribution by Category
sns.kdeplot(data=cleaned_df, x="Contamination", hue="Category", fill=True, alpha=0.5, ax=ax)

# Formatting the plot
#ax.set_title("Distribution of Completeness by Category", fontsize=12, fontweight='bold')
ax.set_xlabel("Contamination (%)", fontsize=10, color='black')
ax.set_ylabel("Density", fontsize=10)
ax.set_ylabel("Density", fontsize=10, color='black')
#ax.legend(title="Category", fontsize=8, title_fontsize=9)
#sns.despine()

# Remove grid lines
ax.grid(False)

# Make the borders black
for spine in ax.spines.values():
    spine.set_edgecolor("black")
    spine.set_linewidth(1.5)

# Save the figure at 600 DPI
save_path = "your_path/Env_MAGS_contamination_.png"    ##### Figure 2E
plt.savefig(save_path, dpi=600, bbox_inches="tight")

# Show the plot
plt.show()

#######################################################################

import matplotlib.pyplot as plt
import seaborn as sns

# Set style without grid
sns.set_style("white")

# Create a figure with an optimized size and high resolution
fig, ax = plt.subplots(figsize=(8, 5))

# Violin plot for Strain Heterogeneity by Category
sns.violinplot(
    data=combined_df,
    x="Category",
    y="Strain_heterogeneity",
    inner="point",
    scale="width",
    palette="viridis",
    linewidth=1.8,  # Make lines sharper
    ax=ax
)

# Formatting the plot
#ax.set_title("Strain Heterogeneity Distribution by Category", fontsize=12, fontweight='bold', color='black')
ax.set_xlabel("", fontsize=10, color='black')
ax.set_ylabel("Strain Heterogeneity (%)", fontsize=10, color='black')

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Remove grid lines
ax.grid(False)

# Make the borders black and slightly thicker
for spine in ax.spines.values():
    spine.set_edgecolor("black")
    spine.set_linewidth(1.5)

# Save the figure at 600 DPI
save_path = "your_path/Env_strain_heterogeneity.png"     ##### Figure 2F
plt.savefig(save_path, dpi=600, bbox_inches="tight")

# Show the plot
plt.show()


##################################################################
import numpy as np

# Set plot style
sns.set_style("whitegrid")

# Create a figure with an optimized size
fig, ax = plt.subplots(figsize=(7, 5))

# Scatter plot for Strain Heterogeneity vs. Contamination
sns.scatterplot(
    data=combined_df,
    x="Contamination",
    y="Strain_heterogeneity",
    hue="Category",
    alpha=0.7,
    edgecolor="black",
    linewidth=0.5,
    ax=ax
)

# Formatting the plot
#ax.set_title("Comparison of Strain Heterogeneity vs. Contamination", fontsize=12, fontweight='bold', color='black')
ax.set_xlabel("Contamination (%)", fontsize=10, color='black')
ax.set_ylabel("Strain Heterogeneity (%)", fontsize=10, color='black')

# Remove grid lines
ax.grid(False)

# Make the borders black and slightly thicker
for spine in ax.spines.values():
    spine.set_edgecolor("black")
    spine.set_linewidth(1.5)

# Save the figure at 600 DPI
scatter_save_path = "your_path/strain_vs_contamination_scatter.png"     ##### Supplementary Figure 1D
plt.savefig(scatter_save_path, dpi=600, bbox_inches="tight")

# Show the scatter plot
plt.show()

########################################
# Compute correlation matrix
correlation_matrix = combined_df[["Contamination", "Strain_heterogeneity"]].corr()

# Plot correlation heatmap
fig, ax = plt.subplots(figsize=(5, 4))
sns.heatmap(correlation_matrix, annot=True, cmap="coolwarm", linewidths=0.5, fmt=".2f", ax=ax)

# Formatting heatmap
#ax.set_title("Correlation Between Strain Heterogeneity & Contamination", fontsize=12, fontweight='bold', color='black')

# Save heatmap at 600 DPI
heatmap_save_path = "your_path/strain_vs_contamination_heatmap.png"     ##### Supplementary Figure 1E
plt.savefig(heatmap_save_path, dpi=600, bbox_inches="tight")  

# Show heatmap
plt.show()

#####################################################

# Import necessary libraries
import matplotlib.pyplot as plt
import seaborn as sns

# Set style (no grid for a cleaner look)
sns.set_style("white")

# Compute category-wise marker lineage counts
marker_category_counts = combined_df.groupby(["Category", "Marker_Lineage"]).size().unstack().fillna(0)

# Get the top 10 most common marker lineages overall
top_marker_lineages = combined_df["Marker_Lineage"].value_counts().head(10).index
filtered_marker_counts = marker_category_counts[top_marker_lineages]

# Create a figure with an optimized size and high resolution
fig, ax = plt.subplots(figsize=(8, 5))

# Plot stacked bar chart
filtered_marker_counts.plot(kind="bar", stacked=True, ax=ax, width=0.7, edgecolor="black")

# Formatting the plot
#ax.set_title("Stacked Bar Chart of Top 10 Marker Lineages by Category", fontsize=12, fontweight='bold', color='black')
ax.set_xlabel("", fontsize=10, fontweight='bold', color='black')
ax.set_ylabel("Count", fontsize=10, color='black')

# Rotate x-axis labels for better readability
plt.xticks(rotation=45, fontsize=9)

# Customize legend (moving it outside for clarity)
legend = ax.legend(title="Marker Lineage", fontsize=7, title_fontsize=8, bbox_to_anchor=(0.72, 1), loc="upper left")
#legend.get_frame().set_edgecolor("black")  # Make legend border black

# Remove grid lines for a cleaner look
ax.grid(False)

# Make the borders black and slightly thicker
for spine in ax.spines.values():
    spine.set_edgecolor("black")
    spine.set_linewidth(1.5)

# Save the figure at 600 DPI
save_path = "your_path/top10_marker_lineages.png"              ### Figure 2G
plt.savefig(save_path, dpi=600, bbox_inches="tight")

# Show the plot
plt.show()



############################################################
import pandas as pd

# Load dataset paths
before_files = {
    "Ambulance": "your_path/PRJNA369713_Fastqscreen_before.csv",
    "Hospital_environment": "your_path/hosp_env_Fastqscreen_before.csv",
    "Hospital_sewage": "your_path/hosp_sewage_Fastqscreen_before.csv",
    "Public_transport": "your_path/public_transport_Fastqscreen_before.csv",
}

after_files = {
    "Ambulance": "your_path/PRJNA369713_Fastqscreen_After.csv",
    "Hospital_environment": "your_path/hosp_env_Fastqscreen_After.csv",
    "Hospital_sewage": "your_path/Hosp_sewage_Fastqscreen_After.csv",
    "Public_transport": "your_path/public_transport_Fastqscreen_After.csv",
}

# Read CSV files into dataframes
before_dfs = {name: pd.read_csv(path) for name, path in before_files.items()}
after_dfs = {name: pd.read_csv(path) for name, path in after_files.items()}

# Define the correct approach: use Read_Processed from each dataset separately
comparison_results = []

for category in before_dfs.keys():
    # Extract before and after data
    before_df = before_dfs[category].copy()
    after_df = after_dfs[category].copy()

    # Ensure column names are correctly formatted (strip spaces)
    before_df.columns = before_df.columns.str.strip()
    after_df.columns = after_df.columns.str.strip()

    # Define relevant metrics for computation
    metrics = ["Unmapped_Human", "Unmapped_Mouse", "Unmapped_Adapters", "Unmapped_Vectors", "Unmapped_Plasmodium"]

    # Compute mapped reads using the correct Read_Processed for each dataset separately
    before_df["Mapped_Human"] = before_df["Read_Processed"] - before_df["Unmapped_Human"]
    before_df["Mapped_Mouse"] = before_df["Read_Processed"] - before_df["Unmapped_Mouse"]
    before_df["Mapped_Adapters"] = before_df["Read_Processed"] - before_df["Unmapped_Adapters"]
    before_df["Mapped_Vectors"] = before_df["Read_Processed"] - before_df["Unmapped_Vectors"]
    before_df["Mapped_Plasmodium"] = before_df["Read_Processed"] - before_df["Unmapped_Plasmodium"]

    after_df["Mapped_Human"] = after_df["Read_Processed"] - after_df["Unmapped_Human"]
    after_df["Mapped_Mouse"] = after_df["Read_Processed"] - after_df["Unmapped_Mouse"]
    after_df["Mapped_Adapters"] = after_df["Read_Processed"] - after_df["Unmapped_Adapters"]
    after_df["Mapped_Vectors"] = after_df["Read_Processed"] - after_df["Unmapped_Vectors"]
    after_df["Mapped_Plasmodium"] = after_df["Read_Processed"] - after_df["Unmapped_Plasmodium"]

    mapped_metrics = ["Mapped_Human", "Mapped_Mouse", "Mapped_Adapters", "Mapped_Vectors", "Mapped_Plasmodium"]

    # Aggregate total mapped reads
    before_totals = before_df[mapped_metrics].sum()
    after_totals = after_df[mapped_metrics].sum()

    # Compute changes in mapped reads
    change = after_totals - before_totals
    percent_change = ((after_totals - before_totals) / before_totals) * 100

    # Create a dataframe for visualization
    category_comparison = pd.DataFrame({
        "Metric": before_totals.index,
        "Before (Mapped)": before_totals.values,
        "After (Mapped)": after_totals.values,
        "Change": change.values,
        "% Change": percent_change.values
    })

    category_comparison["Category"] = category
    comparison_results.append(category_comparison)

# Combine all categories into one dataframe
final_corrected_mapped_comparison_df = pd.concat(comparison_results, ignore_index=True)

# Save the final corrected mapped reads comparison table as CSV
final_corrected_mapped_comparison_df.to_csv("final_corrected_mapped_reads_comparison.csv", index=False)      ### Suppelementary Table 2 (Fast_Screening_result(before and after))

final_corrected_mapped_comparison_df

###########################################################

import matplotlib.pyplot as plt
import seaborn as sns

# Set optimized DPI for high-quality output
optimized_dpi = 600

# Define colors for Before & After
colors = ["#1f77b4", "#ff7f0e"]  # Blue for Before, Orange for After

# Set style for aesthetics
sns.set_style("white")

# Bar Plot for Mapped Reads Before vs. After**
fig, ax = plt.subplots(figsize=(3, 3))

# Extract values for plotting
categories = final_corrected_mapped_comparison_df["Category"].unique()
before_values = final_corrected_mapped_comparison_df.groupby("Category")["Before (Mapped)"].sum()
after_values = final_corrected_mapped_comparison_df.groupby("Category")["After (Mapped)"].sum()

# Define bar positions
x = range(len(categories))
width = 0.5  # Adjust bar width for better spacing

# Plot bars with black borders
ax.bar(x, before_values, width=width, color=colors[0], edgecolor="black", label="Before")
ax.bar([pos + width for pos in x], after_values, width=width, color=colors[1], edgecolor="black", label="After")

# Formatting
#ax.set_title("Mapped Reads Before vs. After Filtering", fontsize=12, fontweight="bold", color="black")
ax.set_xlabel("", fontsize=10, fontweight="bold", color="black")
ax.set_ylabel("Contaminant Reads Count", fontsize=10, color="black")
ax.set_xticks([pos + width /150 for pos in x])
ax.set_xticklabels(categories, rotation=90, fontsize=10, color="black")

# Remove grid lines
ax.grid(False)

# Make the borders black
for spine in ax.spines.values():
    spine.set_edgecolor("black")
    spine.set_linewidth(1.5)

# Add legend
ax.legend(title="", fontsize=10, loc="upper left", frameon=True, edgecolor="black")

# Save the plot
mapped_reads_improved_path = "your_path/Contaminants_reads.png"           #### Figure 2A
plt.savefig(mapped_reads_improved_path, dpi=optimized_dpi, bbox_inches="tight")

# Show the improved plot
plt.show()

#############################################

# Import Required Libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm

# Load Your QUAST Files (Replace with your local file paths)
files = [
    'Ambulance_quast_merged',
    'hosp_env_quast',
    'hosp_sewage_Quast',
    'Public_transport_quast_merged'
]

dfs = [pd.read_csv(f) for f in files]
df = pd.concat(dfs, ignore_index=True)
df.columns = df.columns.str.strip()  # Clean column names

# Convert Key Columns to Numeric
numeric_cols = [
    'No. of contigs', 'Largest contig', 'Total length', 'GC (%)',
    'N50', 'N75', 'L50', 'L75', "No. of N's per 100 kbp"
]
for col in numeric_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')

df = df.dropna(subset=numeric_cols)  # Remove rows with missing values



#  Regression Analysis (Predicting N50)
X = df[['No. of contigs', 'L50', 'Total length', 'Largest contig', 'GC (%)']]
y = df['N50']
X = sm.add_constant(X)
model = sm.OLS(y, X).fit()
print(model.summary())

#  Residual Analysis
df['Predicted N50'] = model.predict(X)
df['Residuals'] = model.resid

# Save regression summary to a text file
regression_report = model.summary().as_text()

# Save to file
report_path = "your_path/Quast_regression_analysis_report.txt"    #### Supplementary Figure 1A
with open(report_path, "w") as file:
    file.write(regression_report)

###################################3###


import matplotlib.pyplot as plt
import seaborn as sns

# Make sure your DataFrame has these columns:
# 'Predicted N50', 'Residuals', 'N50', 'Category'

# Create a 1x2 grid of plots
fig, axes = plt.subplots(1, 2, figsize=(8, 4))  # Side-by-side plots

# Plot 1: Residuals vs. Predicted N50
sns.scatterplot(data=df, x='Predicted N50', y='Residuals', hue='Category', alpha=0.6, ax=axes[0])
axes[0].axhline(0, color='black', linestyle='--')
axes[0].set_title('Residuals vs. Predicted N50')
axes[0].set_xlabel('Predicted N50')
axes[0].set_ylabel('Residuals')
axes[0].legend(title="", fontsize=7, loc="upper left", frameon=True)


# Plot 2: Actual vs. Predicted N50
sns.scatterplot(data=df, x='N50', y='Predicted N50', hue='Category', alpha=0.6, ax=axes[1])
axes[1].plot([df['N50'].min(), df['N50'].max()],
             [df['N50'].min(), df['N50'].max()],
             color='black', linestyle='--', label='Ideal Fit')
axes[1].set_title('Actual vs. Predicted N50')
axes[1].set_xlabel('Actual N50')
axes[1].set_ylabel('Predicted N50')
axes[1].legend(title="", fontsize=7, loc="upper left", frameon=True)

# Final touches
plt.tight_layout()
plt.savefig("Quast_regression_plots.png", dpi=600)  # Supplementary Figure 1C
plt.show()

####################################################

import matplotlib.pyplot as plt
import seaborn as sns

# Your metrics of interest
metrics_to_plot = ['N50', 'No. of contigs', 'L50', 'Total length']

# Set up grid for 2 rows × 3 columns (adjust to fit nicely)
fig, axes = plt.subplots(2, 2, figsize=(6, 8))  # 2 rows, 3 cols
axes = axes.flatten()  # Flatten for easy indexing

# Plot each metric into a subplot
for i, metric in enumerate(metrics_to_plot):
    sns.boxplot(data=df, x='Category', y=metric, ax=axes[i])
    #axes[i].set_title(f'{metric} by Sample Category')
    axes[i].tick_params(axis='x', rotation=90)
    axes[i].set_xlabel('')

# Remove any unused subplot (if we have less than 6)
for j in range(len(metrics_to_plot), len(axes)):
    fig.delaxes(axes[j])

# Final layout and save
plt.tight_layout()
plt.savefig("Quast_combined_boxplots.png", dpi=600)    # Supplementary Figure 1B
plt.show()


#########################################################

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats

# Set figure size
plt.figure(figsize=(3.3, 4))

# Apply log scale to X-axis for better visualization
plt.xscale("log")

# Scatter plot with transparent, smaller points
sns.scatterplot(data=df, x='No. of contigs', y='N50', hue='Category', alpha=0.6, s=30)

# Fit and plot a single regression line
sns.regplot(data=df, x='No. of contigs', y='N50', scatter=False, color='black', logx=True, line_kws={'linewidth': 3})

# Compute and display R²
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(df['No. of contigs']), df['N50'])
plt.text(0.96, 0.93, f'R² = {r_value**2:.2f}', transform=plt.gca().transAxes,
         fontsize=8, ha='right', bbox=dict(facecolor='white', alpha=0.1))

# Titles and labels
plt.xlabel('Number of Contigs', fontsize=10)
plt.ylabel('N50', fontsize=10)

# Improve legend placement
plt.legend(title="", fontsize=6, frameon=True, loc='center',  bbox_to_anchor=(0.74, 0.80))

# **Increase Border Thickness**
for spine in plt.gca().spines.values():
    spine.set_linewidth(1.5)  # Set all borders to thickness 1.5

# Save and show
plt.tight_layout()
plt.savefig("n50_vs_contigs_final.png", dpi=600)       #### Figure 2B
plt.show()

