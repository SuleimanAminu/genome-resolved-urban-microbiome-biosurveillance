# -----------------------------------------------
# Script: 03_subset_pathogenic_species.py
# Authors: Suleiman & AbdulAziz
# Purpose: Extract pathogenic species (CARD + VFDB) from the corrected abundance matrix
# -----------------------------------------------

import pandas as pd

# Load the manually curated list of pathogenic species from CARD + VFDB
pathogenic_species = pd.read_csv('/content/vfdb_card_matched_species')

# Load the corrected species abundance matrix
species_final = pd.read_csv('/content/species_final.csv')   ###### make sure the first column is named "taxon" before running

# --- Data Cleaning ---
# Rename first column in the pathogen list to 'taxon' for clarity
pathogenic_species = pathogenic_species.rename(columns={pathogenic_species.columns[0]: "taxon"})

# Strip whitespace for robust matching
pathogenic_species["taxon"] = pathogenic_species["taxon"].str.strip()
species_final["taxon"] = species_final["taxon"].str.strip()

# --- Subset species_final to only include pathogenic species ---
species_list = pathogenic_species["taxon"].tolist()
filtered_data = species_final[species_final["taxon"].isin(species_list)]

# --- Save the filtered pathogenic species matrix ---
output_path = '/content/Env_Pathogenic_species.csv'
filtered_data.to_csv(output_path, index=False)

print(f"Pathogenic species matrix saved to: {output_path}")
print(filtered_data.head())

