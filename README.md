
#  Genome-Resolved Urban Microbiome Biosurveillance

A fully integrated framework for **genome-resolved metagenomics**, **ecological modeling**, and **machine learning** designed to support microbial surveillance and risk assessment in urban environments.

This pipeline enables high-resolution tracking of pathogen-associated taxa across infrastructure types such as ambulances, hospital interiors, sewage systems, and public transport. It combines robust preprocessing, species-level ecological analytics, and predictive modeling to uncover microbial transmission patterns, contamination risks, and ecological vulnerability.

---

## What This Framework Offers

- **Genome-resolved profiling**: From raw reads to annotated MAGs with virulence and AMR gene screening.
-  **Ecological modeling**: Diversity analysis, indicator species detection, community structuring.
- **Machine learning diagnostics**: Environmental source prediction, entropy-based identity modeling, and composite risk scoring.
- **Simulation-based inference**: Synthetic community blending, ecological fragility assessment, and uncertainty quantification.

Built for reproducibility, explainability, and cross-environment comparability, this framework operationalizes biosurveillance as a diagnostic, data-driven process.

---


## Modules

### 1 Bioinformatics
High-throughput metagenomic processing on HPC:
- Quality control, adapter trimming
- Host decontamination and screening
- Assembly, binning, MAG quality control
- Functional annotation (VFDB & CARD)

**Scripts:**
```
01_Bioinformatics/
├── 01_run.sh
├── 02_fastq_screen.sh
├── 03_collect_files.sh
├── simplifyFastaHeaders.pl
```

---

### 2 Quality & Batch Subsetting
Species-level matrix refinement and pathogen-focused subsetting:
- MAG quality metrics & contiguity modeling
- Batch correction (PCA, UMAP, t-SNE)
- Pathogen selection based on virulence and resistance

**Scripts:**
```
02_Quality_batch_subsetting/
├── 04_mag_quality_metrics_analysis.py
├── 05_normalize_species_counts.py
├── 06_Batch_Correction_PCA_UMAP_tsne.R.R
├── 07_subset_pathogenic_species.py
```

---

### 3 Ecology & Diversity Analysis
Ecological fingerprinting of urban environments:
- Alpha & beta diversity, PERMANOVA, ANOSIM
- Indicator species and NMDS ordination
- WHO-pathogen co-occurrence networks

**Scripts:**
```
03_Ecology/
├── 08_species_community_analysis.R
├── 09_pathogen_prevalence_diversity.py
```

---

### 4 Machine Learning & Risk Modeling
Explainable AI pipeline for microbial classification:
- Random Forest classifier training and benchmarking
- Monte Carlo simulation with feature importance
- Synthetic community blending
- Entropy modeling and risk scoring

**Scripts:**
```
04_Machine_Learning/
├── 10_rf_model_training.py
├── 11_model_comparison.py
├── 12_ML_Comparison_Script.R
├── 13_ML_Env_Biosurveillance.py
```

---

## Requirements

- Python ≥ 3.8
- R ≥ 4.3.0
- Conda (recommended)
- SLURM-compatible HPC environment (for `.sh` scripts)

Install Python environment:
```bash
conda create -n urban_env python=3.11
conda activate urban_env
pip install pandas numpy scikit-learn seaborn matplotlib joblib tqdm
```

---

## Usage Examples

Run ML classifier:
```bash
sbatch 04_Machine_Learning/10_rf_model_training.py
```

Run ecology analysis:
```bash
Rscript 03_Ecology/08_species_community_analysis.R
```

---

##  Folder Structure

```
Urban_Microbiome_Surveillance/
├── 01_Bioinformatics/
├── 02_Quality_batch_subsetting/
├── 03_Ecology/
├── 04_Machine_Learning/
├── README.md
├── LICENSE
└── .gitignore
```

---

## License

MIT License — free to use, adapt, and cite with attribution.

---

## Citation

If you use this Framework, please cite:
**Aminu S.**, Ascandari A., Benhida R., Daoud R. (2025).  
*GRUMB: Genome-Resolved Urban Microbiome Biosurveillance*.  
Zenodo. [https://doi.org/10.5281/zenodo.15505402](https://doi.org/10.5281/zenodo.15505402)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15505402.svg)](https://doi.org/10.5281/zenodo.15505402)



## Submitted Articles Related to the Framework

Aminu, S., Ascandari, A., Mokhtar, M.M., El Allali, A., Benhida, R., Daoud, R. (2025).  
*Genome-Resolved Species-Level Surveillance and Predictive Risk Modeling of Urban Microbiomes*.  
Microbiome** (Submitted).

Aminu, S., Ascandari, A., Benhida, R., Daoud, R. (2025).  
*GRUMB: A Genome-Resolved Metagenomic Framework for Monitoring Urban Microbiome and Diagnosing Pathogen
Risk*.  
**Bioinformatics** (Submitted).
 

