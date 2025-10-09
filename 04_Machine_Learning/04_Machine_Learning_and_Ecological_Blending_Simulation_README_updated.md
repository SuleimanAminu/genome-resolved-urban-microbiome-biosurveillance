

# *04_Machine_Learning_and_Ecological_Blending_Simulation Module*

This module implements **machine learning–based ecological inference** for urban microbiome datasets, providing predictive modeling, feature stability analysis, and ecological simulation frameworks.  
Together, these scripts quantify **classification robustness**, **feature-level ecological signatures**, and **synthetic contamination resilience** across environmental microbiomes.

---

## **Scripts**
```
04_Machine_Learning_and_Ecological_Blending_Simulation/
├── 01_Model_selection_Holdout_Nested.py
├── 02_Feature_importance_stability.py
├── 03_Contaminations_Simulations_Ecological_Scoring.py

```
---

## **Scripts and Functions** 

## `01_Model_selection_Holdout_Nested.py` — Model Benchmarking & Nested Cross-Validation

## **Overview**
Performs **supervised machine learning classification** using microbial abundance data to discriminate between urban environments (e.g., hospital, sewage, ambulance, public transport).  
This script evaluates and statistically compares multiple algorithms using **cross-validation, nested tuning, and hold-out testing**.

## **Key Analyses**
- **Model Benchmarking:**
  - Classifiers: Random Forest, Decision Tree, SVM, Logistic Regression, Gradient Boosting  
  - Metrics: Accuracy, Balanced Accuracy, F1-macro, Cohen’s Kappa  
  - 10-fold stratified CV with performance mean ± SD  
- **Model Statistics:**
  - Friedman test, Wilcoxon posthoc comparisons, and Nemenyi heatmaps for model significance  
- **Nested CV + Holdout:**
  - 5×3 nested cross-validation for Decision Tree & Random Forest  
  - Saves fold-level metrics and final tuned hyperparameters  
  - Evaluates final model on hold-out data (20%) with confusion matrices and per-class F1/Kappa  
- **Visualizations:**
  - Model benchmarking plots with error bars  
  - Confusion matrices and NestedCV vs Hold-out accuracy comparisons  

## **Outputs**
- `ML_Benchmarking_AllModels_1.csv`, `ML_FriedmanTest.csv`  
- `NestedCV_summary.csv`, `Holdout_metrics.csv`, `Holdout_confusion_matrix.csv`  
- Figures:
  - `ML_Benchmarking_1.tiff`, `Figure_ConfusionMatrix_Counts.tiff`  
  - `Figure_NestedCV_vs_Holdout.tiff`, `ML_NemenyiHeatmap_*.tiff`

## **Scientific Context**
Establishes the **best-performing classifier** for environment prediction and ensures robust generalization via nested CV, forming the foundation for feature importance and contamination simulation analyses.

---

## `02_Feature_importance_stability.py` — Feature Importance & Stability Analysis

## **Overview**
Quantifies the **stability, reproducibility, and consensus ranking** of microbial features driving classification accuracy across repeated Random Forest and Decision Tree models.

## **Key Analyses**
- **Feature Importance:**
  - Repeated stratified 5-fold CV (50 iterations)  
  - Random Forest Gini and permutation importance for top-K features (default = 60)  
  - Partial saves (`importance_partials/`) for reproducibility  
- **Consensus Ranking:**
  - Integrates importance scores across RF, DT, and permutation methods  
  - Produces normalized consensus feature ranking (`Consensus_Feature_Ranking.csv`)  
- **Stability Visualization:**
  - Boxplots for top-15 species (importance distribution)  
  - Highlighted barplot of top-30 consensus features  
- **Per-Environment Models:**
  - One-vs-rest RF classifiers to identify **environment-specific top species**  
  - Generates per-environment feature barplots and top species tables  
- **Ecological Projection:**
  - Heatmap of top 50 discriminative species (z-scored across environments)  

## **Outputs**
- `RF_importances_repeated.csv`, `Permutation_importances_topK.csv`  
- `Consensus_Feature_Ranking.csv`  
- Figures:
  - `Consensus_Top30_Barplot.tiff`, `RF_ImportanceStability_Top15.tiff`  
  - `TopFeatures_EnvSpecific_Barplots.tiff`, `TopSpecies_EnvHeatmap.tiff`

## **Scientific Context**
Provides a **feature-level ecological fingerprint**, pinpointing microbial taxa that consistently define environmental identity and robustness.

---

## `03_Contaminations_Simulations_Ecological_Scoring.py` — Ecological Simulations & Risk Modeling

## **Overview**
Implements the **Synthetic Mixing Framework** to simulate contamination events between environments, quantify **Minimal Detectable Contamination (MDC)**, and evaluate **ecological resilience and risk** using entropy and dominance metrics.

## **Key Analyses**
1. **Dose–Response Simulation:**
   - Linear blending of donor–recipient microbial profiles (`α = 0–1`)  
   - Computes probability of donor classification and confidence intervals  
   - Derives MDC values (threshold α for ≥80% donor prediction)  
2. **Entropy-Based Stability:**
   - Shannon entropy across class probabilities  
   - Identifies **tipping points** of ecological identity collapse  
   - Switch-point extraction for directional dominance analysis  
3. **Dominance & Network Matrices:**
   - Heatmaps showing cross-environment dominance (wins/losses)  
   - Captures compositional hierarchy and directional influence  
4. **Ecological Risk Scoring:**
   - Integrates metrics: Overridden frequency, Entropy, Richness, Variance  
   - Computes normalized risk under three weighting schemes:  
     - Weighted primary (0.4–0.3–0.15–0.15)  
     - Equal weights  
     - PCA-derived weights  
   - Radar and bar plots visualize comparative environmental risks  

## **Outputs**
- `Simulations_donor_recipient_withCI.csv`, `Simulations_MDC_table_refined.csv`  
- `synthetic_blend_entropy.csv`, `synthetic_transition_points.csv`  
- `environment_risk_scores_all.csv`, `risk_score_rank_correlations.csv`  
- Figures:
  - `Figure_Sim_MDC_heatmap_sample.tiff`, `EntropyCurve_*.tiff`  
  - `Environment_Dominance_Matrix.tiff`, `Contamination_Risk_Barplot.tiff`  
  - `Ecological_Risk_Radar_Clean.tiff`

## **Scientific Context**
Links machine learning predictions with **ecological theory**, providing quantifiable measures of:
- Microbial identity resilience  
- Contamination susceptibility  
- Directional dominance  
- Ecological risk hierarchies  

---

##  **Output Summary** 

| Script | Focus | Outputs |
|--------|--------|----------|
| `01_Model_selection_Holdout_Nested.py` | Model benchmarking & validation | CV results, holdout metrics, confusion matrices |
| `02_Feature_importance_stability.py` | Feature importance reproducibility | Consensus ranks, per-environment top species |
| `03_Contaminations_Simulations_Ecological_Scoring.py` | Synthetic mixing & ecological risk | Dose–response, entropy, dominance, and risk plots |

---

##  **Requirements** 

- **Python (≥ 3.8)**  
  Packages:  
  `pandas`, `numpy`, `scikit-learn`, `scipy`, `matplotlib`, `seaborn`, `scikit-posthocs`, `Pillow`

---

## **Example Usage** 

```bash
# 1. Model benchmarking and nested validation
python 01_Model_selection_Holdout_Nested.py

# 2. Feature stability and consensus ranking
python 02_Feature_importance_stability.py

# 3. Synthetic contamination simulations and ecological scoring
python 03_Contaminations_Simulations_Ecological_Scoring.py

```
