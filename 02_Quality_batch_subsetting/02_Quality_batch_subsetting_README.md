
# Quality Correction & Pathogen Subsetting Module

This module provides post-assembly normalization, batch correction, and selective filtering of pathogenic species from species-level metagenomic abundance data. It acts as a critical preprocessing stage before ecological analysis or machine learning.

---

##  Scripts and Functions

### `04_mag_quality_metrics_analysis.py` — MAG Quality Metrics & Technical Modeling

#### Overview:
Analyzes MAG metrics across environmental samples to assess:
- Completeness, contamination, strain heterogeneity
- QUAST assembly statistics (N50, L50, total length, GC%)
- Marker lineage diversity and comparison across environments
- Linear regression to model N50 using technical variables

 **Note:** Ensure required columns (e.g., `Completeness`, `Contamination`, `Category`) are present in your MAG summary tables.

---

### `05_normalize_species_counts.py` — Taxonomic Matrix Normalization

#### Overview:
Converts Bracken outputs into a unified, scaled species abundance matrix.

#### Key Steps:
- Reads raw Bracken files and scales by read depth
- Normalizes per-sample abundances
- Merges all environments into one matrix (`Merged_species.csv`)
- Generates final pivoted matrix (`Species_expression_matrix.csv`)
- Aligns with curated metadata (`Species_metadata.csv`)


---

### `06_Batch_Correction_PCA_UMAP_tsne.R.R` — Dimensionality Reduction & Correction

#### Overview:
Removes batch effects and visualizes community clustering using:
- PCA
- UMAP
- t-SNE
- PERMANOVA diagnostics

#### Pipeline:
1. Load matrix + metadata
2. Hellinger transformation (`vegan::decostand`)
3. PCA, UMAP, and t-SNE visualization
4. Batch correction with `limma::removeBatchEffect`
5. PERMANOVA before/after correction
6. Export corrected matrix and plots

 **Requires R packages:** `vegan`, `ggplot2`, `lme4`, `limma`, `umap`, `Rtsne`, `shadowtext`, etc.

---

### `07_subset_pathogenic_species.py` — Pathogen-Associated Species Extraction

#### Overview:
Filters `species_final.csv` using a curated list of species carrying:
- ARGs (from CARD)
- Virulence factors (from VFDB)

#### Output:
- `Env_Pathogenic_species.csv`: Contains only clinically relevant species.

Useful for:
- Machine learning on pathogenic fingerprints
- Entropy-based identity fragility scoring
- AMR surveillance and risk prioritization

---


## Requirements

- Python ≥ 3.7
- pandas ≥ 1.0
- R ≥ 4.0
- Required R packages listed above

---

## Example Usage

```bash
python 05_normalize_species_counts.py
Rscript 06_Batch_Correction_PCA_UMAP_tsne.R
python 07_subset_pathogenic_species.py
```

---

##  Scientific Relevance

This module prepares clean, interpretable species-level data for:
- Ecological ordination (e.g., NMDS, PERMANOVA)
- Machine learning classifiers
- Cross-environmental comparisons
- Targeted biosurveillance and risk profiling
