
# Machine Learning & Risk Modeling Module

This module implements supervised classification, feature importance analysis, and entropy-based ecological risk modeling using species-level abundance data from pathogenic microbes in urban environments.

---

##  Scripts and Functions

### `10_rf_model_training.py` — Random Forest Training with GridSearchCV

#### Overview:
Trains a Random Forest classifier to predict the environmental origin (e.g., ambulance, sewage) of a sample based on pathogenic species profiles.

#### Key Features:
- Uses `GridSearchCV` for hyperparameter optimization
- Saves trained model and scoring metrics
- Supports scalable training on HPC clusters

Submit via SLURM using `run_rf_training.sh`:
```bash
sbatch run_rf_training.sh
```

---

### `model_comparison.py` — Multi-Classifier Benchmarking

#### Overview:
Evaluates six machine learning models using 10-fold stratified cross-validation.

#### Metrics:
- Accuracy
- Macro F1-score
- Cohen’s Kappa

Supports identification of the most generalizable model.

Submit via SLURM using `model_comparison.sh`:
```bash
sbatch model_comparison.sh
```

 Visualization of comparison results can be done using:
- `ML_Comparison_Script.R` (R script)

---

### `13_ML_Env_Biosurveillance.py` — Monte Carlo Learning & Risk Simulation

#### Overview:
Implements a complete ML ecosystem for ecological diagnostics:
- Monte Carlo learning for performance and feature stability
- Synthetic microbial mixing simulations
- Entropy and override-based risk scoring

#### Workflow Steps:

**1. Data Preparation**
- Joins abundance matrix with metadata
- Reshapes to long-form
- Drops technical metadata (e.g., Instrument)

**2. Monte Carlo Simulation**
- Trains 100 Random Forest models with random 80/20 splits
- Records prediction accuracy, F1, Kappa, confusion matrices

**3. Feature Importance**
- Extracts per-species importance scores across runs
- Outputs mean and variance (`monte_carlo_species_importance.csv`)

**4. Synthetic Microbial Blending**
- Blends species matrices from all environment pairs at 10% intervals
- Tracks prediction switches, confidence scores, and entropy

**5. Ecological Risk Scoring**
- Computes:
  - Prediction override frequency
  - Classification uncertainty (Shannon entropy)
  - Intra-environmental variance
  - Species richness
- Outputs final `Contamination_Risk_Score.csv`

---

##  Output Files

| Filename                               | Description                                   |
|----------------------------------------|-----------------------------------------------|
| `trained_model.pkl`                    | Final optimized RF model                      |
| `mc_confusion_matrices.csv`           | Confusion matrices for each simulation        |
| `monte_carlo_species_importance.csv`   | Feature importance across 100 RF models       |
| `entropy_scores.csv`                   | Entropy by blending ratio and environment     |
| `Contamination_Risk_Score.csv`         | Composite ecological vulnerability score      |

---

## Requirements

- Python ≥ 3.8
- Required packages: `scikit-learn`, `pandas`, `numpy`, `joblib`, `matplotlib`, `tqdm`, `seaborn`

---

##  Example Usage

```bash
sbatch run_rf_training.sh
sbatch model_comparison.sh
python 13_ML_Env_Biosurveillance.py
```

---

## Scientific Relevance

This module bridges:
- Predictive classification of microbiomes
- Explainable feature attribution for sentinel species
- Simulation-based inference of ecological vulnerability

Supports scalable microbial risk stratification in line with WHO and One Health surveillance priorities.
