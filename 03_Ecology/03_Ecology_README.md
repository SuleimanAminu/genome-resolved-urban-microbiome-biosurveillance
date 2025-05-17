
# Ecology & Community Structure Module

This module performs comprehensive ecological analysis on batch-corrected species-level microbial abundance data, focusing on diversity, indicator species, and the ecological behavior of WHO-priority pathogens in urban infrastructure environments.

---

## Scripts and Functions

### `08_species_community_analysis.R` — Diversity, Indicator Species & Network Inference

#### Overview:
This R script provides in-depth ecological analytics, including:
- Alpha diversity metrics: Shannon, Simpson, and species richness
- Kruskal–Wallis tests and Dunn’s post hoc comparisons
- Beta diversity:
  - PERMANOVA (adonis2)
  - NMDS ordination
  - ANOSIM and multivariate dispersion
- Prevalence classification:
  - Core (>80%)
  - Secondary (50–80%)
  - Peripheral (<50%)
- Indicator species analysis (IndVal)
- WHO-pathogen co-occurrence network inference

 Outputs:
- High-resolution TIFF figures for publication
- Supplementary tables with p-values, R², and indicator values
- Ecological fingerprints of each infrastructure environment

️ **Important Note:**
Avoid loading `igraph` until network analysis — it masks `union()`, `diversity()` and other functions used in `vegan` and `dplyr`.

---

### `09_pathogen_prevalence_diversity.py` — WHO Pathogen Prevalence & Impact on Diversity

#### Overview:
Analyzes the occurrence and ecological impact of WHO-priority pathogens across urban microbial communities.

#### Key Functions:
- Prevalence analysis for exact and genus-level WHO-pathogens
- Stacked bar charts of top 30 pathogens by environment
- Grouped violin plots of pathogen distributions
- Shannon diversity comparison:
  - With vs. without each pathogen
  - T-tests for statistical difference
  - Boxplots showing pathogen-specific diversity shifts

#### Code Note:
Diversity is computed using:
```python
abundance_matrix['Shannon_Diversity'] = abundance_matrix.apply(lambda row: entropy(row + 1e-9, base=np.e), axis=1)
```

 This script supports both **epidemiological risk** and **ecological interpretation** of pathogen behavior.

---

## Output Files

- `diversity_summary.csv`: All alpha diversity scores per sample
- `indicator_species.csv`: High-confidence environment markers
- `nmds_coords.csv`, `permanova_summary.csv`: Beta diversity results
- `who_pathogen_prevalence.csv`: Environment × pathogen presence matrix
- `shannon_comparison.csv`: Diversity impact of pathogen presence

---

##  Requirements

- **R (≥ 4.4)**: `vegan`, `ggplot2`, `indicspecies`, `dplyr`, `car`, `shadowtext`
- **Python (≥ 3.8)**: `pandas`, `scipy`, `matplotlib`, `numpy`

---

##  Example Usage

```bash
Rscript 08_species_community_analysis.R
python 09_pathogen_prevalence_diversity.py
```

---

##  Scientific Relevance

This module enables:
- Ecological fingerprinting of infrastructure-specific microbial communities
- Quantitative profiling of high-risk environments
- Identification of sentinel taxa for biosurveillance
- Integration of AMR pathogen dynamics with diversity theory

Suitable for:
- WHO-aligned risk mapping
- Microbial ecosystem comparison
- Epidemiological enrichment studies
