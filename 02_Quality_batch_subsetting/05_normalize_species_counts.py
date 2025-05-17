# Importing Panda
import pandas as pd

###############################################################
# STEP 1: Normalize Bracken species-level data per environment
###############################################################

# NOTE: Repeat this step separately for each environment:
# i.e., Ambulance_species.csv, Public_transport_species.csv, hosp_env_Species.csv

df = pd.read_csv('your_path/hosp_sewage_Species.csv')  # Replace filename as needed

# Total read count per sample (needed for scaling)
df['total_reads_in_sample'] = df.groupby('Sample_ID')['new_est_reads'].transform('sum')

# Compute scaled counts using fraction and total read depth
df['scaled_counts'] = (df['fraction_total_reads'] * df['total_reads_in_sample']).round().astype(int)

# Aggregate data by taxonomy per sample
aggregated_df = df.groupby(
    ['Sample_ID', 'Taxonomy_name', 'Taxonomy_id', 'Taxonomy_lvl', 'Category']
).agg({
    'scaled_counts': 'sum',
    'total_reads_in_sample': 'first'
}).reset_index()

# Add relative abundance
aggregated_df['relative_abundance'] = aggregated_df['scaled_counts'] / aggregated_df['total_reads_in_sample']

# Save corrected output for each environment
aggregated_df.to_csv('your_path/Corrected_hosp_sewage_species.csv', index=False)


###############################################################
# STEP 2: Merge all normalized environments into one file
###############################################################

# Load corrected species tables for all environments
df1 = pd.read_csv('your_path/Corrected_Ambulance_species.csv')
df2 = pd.read_csv('your_path/Corrected_Public_transport_species.csv')
df3 = pd.read_csv('your_path/Corrected_hosp_env_species.csv')
df4 = pd.read_csv('your_path/Corrected_hosp_sewage_species.csv')

# Combine into one long-format dataframe
merged_df = pd.concat([df1, df2, df3, df4], ignore_index=True)

# Save merged species table
merged_df.to_csv('your_path/Merged_species.csv', index=False)
print("Merged species data:", merged_df.shape)


###############################################################
# STEP 3: Generate Expression Matrix & Align Metadata
###############################################################

# Load merged species data and metadata
data = pd.read_csv('your_path/Merged_species.csv')
metadata = pd.read_csv('your_path/Merged_metadata_all.csv')

# Merge on Sample_ID to integrate group/environmental attributes
merged_data = pd.merge(data, metadata, on='Sample_ID')

# Optional: Aggregate counts if duplicate taxonomies exist per sample
data_aggregated = merged_data.groupby(
    ['Taxonomy_name', 'Sample_ID'], as_index=False
)['scaled_counts'].mean()

# Create species expression matrix (Taxa x Samples)
expression_matrix = data_aggregated.pivot(
    index='Taxonomy_name',
    columns='Sample_ID',
    values='scaled_counts'
).fillna(0)

# Prepare metadata matrix (Sample_ID as index)
metadata_relevant = merged_data[[
    'Sample_ID', 'Project_ID', 'Group', 'Country',
    'environmental_material', 'Region', 'Continent',
    'Instrument', 'Sequencing_Center'
]].drop_duplicates().set_index('Sample_ID')

# Ensure matrices are aligned
expression_matrix = expression_matrix[metadata_relevant.index]

# Report shapes for sanity check
print("Expression Matrix Shape:", expression_matrix.shape)
print("Metadata Matrix Shape:", metadata_relevant.shape)

# Save final matrices
expression_matrix.to_csv('your_path/Species_expression_matrix.csv')
metadata_relevant.to_csv('your_path/Species_metadata.csv')

print("Species expression matrix and metadata saved successfully.")

