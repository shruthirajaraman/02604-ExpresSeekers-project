# ExpresSeekers
# Lung and Breast Cancer RNA-Seq Analysis

This repository contains the complete analysis pipeline for identifying differentially expressed genes (DEGs), performing gene ontology (GO) enrichment, and training classification models using RNA-seq data from cancer datasets.

## Folder Overview

### DataPreprocessing  
Scripts for cleaning and formatting raw RNA-seq data.

- `DataPrepForDESeq2.ipynb` / `DataPrepForDESeq2.py`  
  Prepare RNA-seq count data from GDC for DESeq2. Both files perform the same preprocessing steps; choose the notebook or script as needed.

### DataAnalysis  
Contains statistical analyses, including exploratory data analysis, PCA, and silhouette score evaluations.

### DESeq2  
DESeq2-based differential expression analysis of cancer subtypes vs. normal tissue, with code adapted for RNA-seq count data from GDC.

### GOAnalysis  
Performs gene ontology enrichment (biological process, molecular function, KEGG pathways) on DEG sets.

- `GeneOntology.py`  
  Runs GO enrichment analysis and visualizes significant terms using bar plots and colormaps of combined scores across GO terms.

### ModelTraining  
Supervised learning models (e.g., XGBoost, SVM, Logistic Regression) trained on selected gene sets to classify lung cancer subtypes.

## Contributors
*   [Lalithashruthi Rajaraman](https://github.com/shruthirajaraman)
*   [Simran Sodhi](https://github.com/Simran-Sodhi)
*   [Sanchitha Kuthethoor](https://github.com/SanchithaK)
