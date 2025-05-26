# Vertical Integration of CITE‑seq Data

This repository contains code and results for comparing three methods of vertically integrating CITE‑seq data (RNA + ADT): Principal Component Analysis (PCA), Seurat WNN (v5.3.0), and scArches TOTALVI (v0.6.1).

---

## Table of Contents

- [Overview](#overview)  
- [Data](#data)  
- [Methods](#methods)  
  - [1. PCA](#1-pca)  
  - [2. Seurat WNN](#2-seurat-wnn)  
  - [3. scArches TOTALVI](#3-scarches-totalvi)  
- [Benchmark Metrics](#benchmark-metrics)  


---

## Overview

CITE‑seq measures both transcriptome (RNA) and surface protein (ADT) in the same single cells. Vertical integration methods aim to combine these two modalities into a joint latent representation for downstream tasks such as clustering and visualization. Here, we compare:

1. **PCA** on concatenated RNA + ADT features  
2. **Seurat WNN** (Weighted Nearest Neighbors)  
3. **scArches TOTALVI** (Variational Autoencoder)

Each method’s learned latent embeddings are visualized via UMAP and evaluated with quantitative benchmarks.

---

## Data

- **Input**  
  - RNA count matrix (cells × genes)  
  - ADT count matrix (cells × proteins)  

- **Preprocessing**  
  - Normalize and log‑transform RNA  
  - Normalize ADT  
  - Select 4,000 highly variable genes for RNA  

---

## Methods

### 1. PCA

1. Concatenate RNA and ADT feature matrices (cells × (genes + proteins)).  
2. Run PCA to obtain a 50‑dimensional latent space.  
3. Save `pca_latent.csv` and UMAP coordinates `pca_umap.csv`.

### 2. Seurat WNN

1. Process RNA and ADT separately:  
   - Normalize  
   - Identify variable features  
   - Scale and run PCA  
2. Integrate modalities with WNN (Seurat v5.3.0).  
3. Compute UMAP on WNN graph and cluster with Leiden.  
4. Save `wnn_latent.csv` and `wnn_umap.csv`.

### 3. scArches TOTALVI

1. Preprocess RNA as above (normalize, log‑transform, variable genes).  
2. Fit TOTALVI model (scArches v0.6.1) to learn joint latent representation.  
3. Extract 50‑dimensional latent embeddings and UMAP.  
4. Save `totalvi_latent.csv` and `totalvi_umap.csv`.

---

## Benchmark Metrics

We assessed integration quality using **scib** (v1.17) metrics:

| Metric                | Description                                                  |
|-----------------------|--------------------------------------------------------------|
| **ARI**               | Adjusted Rand Index (clustering agreement)                   |
| **NMI**               | Normalized Mutual Information (cluster label agreement)      |
| **cLISI_Graph**       | Cell‑type Local Inverse Simpson’s Index (graph‑based)        |
| **Silhouette**        | Silhouette score on latent embedding                         |
| **Isolated Labels ASW** | Average silhouette width for isolated cell‑type labels     |

Scripts in `benchmark/` compute these metrics and output `benchmark_results.csv`.

---

