# metatest_final — Multi-Omics Integration Framework (MOFA2 + Biomarker Selection)

## Overview

This repository contains the **structured version** of my multi-omics integration framework developed during my PhD research.
It integrates multiple omics layers—**proteomics, metabolomics, transcriptomics and clinical variables**—using latent factor modelling (MOFA2) and a multi-step biomarker selection strategy.

---

## Main Features

### 1. **Multi-omics input layers**

* **Proteins**
* **Metabolites**
* **Transcriptomics**
* **Clinical features / phenotypes**

Preprocessing and harmonization scripts are included to ensure consistent sample matching across all modalities.

---

### 2. **Latent factor modelling with MOFA2**

* Construction of multi-omics models
* View-specific feature selection
* Factor interpretation and visualization
* Hyperparameter exploration
* Projection of new data onto latent spaces

MOFA2 serves as the **core integrative model**, enabling dimensionality reduction across omics and identification of shared and modality-specific variation.

---

### 3. **Biomarker selection pipeline**

The project includes several complementary layers of feature selection, such as:

* **MOFA2 weight thresholds** (per view)
* **sPLS-DA feature stability**
* **LASSO multinomial logistic regression**
* **Multinomial modelling via nnet::multinom**
* **BRMS Bayesian categorical modelling**
* **Ranking by metrics:** Balanced Accuracy, LogLoss, AUC
* **Extraction of biomarker coefficients via `fixef()`**

These steps form a **multi-criteria biomarker discovery framework**, combining unsupervised latent factors with supervised classification models.

---

### 4. **Model evaluation**

Includes:

* Cross-validation loops
* Repeated seeds
* ROC curves (multi-class)
* Feature stability frequency tables
* Final ranking of biomarkers

---
## Rationale

This repository consolidates the **reproducible version** of a multi-omics workflow that evolved through exploratory MSc work (integromics_2) and extensive PhD refinement.
It reflects:

* real multi-omics complexity,
* iterative model improvement,
* integration of modern statistical and Bayesian methods,
* and practical biomarker discovery needs.

---

## Status

The pipeline is functional and has been used for large-scale PCOS multi-omics analyses.
While there is room for further modularization, this version represents the **most stable and interpretable** stage of the project.

