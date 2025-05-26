import scanpy as sc
import anndata as ad
import sys
import torch
import numpy as np
import scarches as sca
import matplotlib.pyplot as plt
import numpy as np
import scvi as scv
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
import gdown
import json
import time
import os
import anndata as ad
from scipy.io import mmwrite

import warnings
warnings.filterwarnings("ignore")

# Define paths
RNA_path = "/public/home/chenjiaminggroup/shenibdproject/multimodal/GSM4972212/"

# Read RNA data
rna_data = pd.read_csv(RNA_path + 'GSM4972212_Citeseq_Human.GBM.R2_5_ND8.filtered.RNA.feature.bc.matrix.csv.gz', index_col=0).T
cell_names = rna_data.index  # Barcodes
gene_names = rna_data.columns  # Features
M_rna = csr_matrix(rna_data.values)

# Read ADT data
adt_data = pd.read_csv(RNA_path + 'GSM4972212_Citeseq_Human.GBM.R2_5_ND8.filtered.ADT.feature.bc.matrix.csv.gz', index_col=0).T
adt_features = adt_data.columns  # ADT features
M_adt = csr_matrix(adt_data.values)

# Read annotation data (if needed for downstream analysis)
annot = pd.read_csv(
    RNA_path + 'GSM4972212_annot.Citeseq_Human.GBM.R2_5_ND8.csv.gz',
    sep=',',
    dtype=str,             # 全部先当作字符串
    quotechar='"'
)
# annot
annot = annot.set_index('cell')
annot = annot.rename(columns={'cluster': 'cell_type'})


# Create AnnData object with RNA data
adata_RNA = ad.AnnData(M_rna, obs=pd.DataFrame(index=cell_names), var=pd.DataFrame(index=gene_names))
adata_RNA.obs['cell_type'] = annot['cell_type']
adata_RNA.obs
adata_RNA.var_names_make_unique()
adata_RNA.layers["counts"] = adata_RNA.X.copy()



# Preprocess RNA data
sc.pp.normalize_total(adata_RNA)
sc.pp.log1p(adata_RNA)
sc.pp.highly_variable_genes(
    adata_RNA,
    flavor="seurat_v3",
    n_top_genes=4000,
    subset=False
)
adata_RNA = adata_RNA[:, adata_RNA.var.highly_variable].copy()
adata_RNA.obsm["protein_expression"] = adt_data


# Add batch information
modality = ['batch1'] * adata_RNA.shape[0]
adata_RNA.obs['batch'] = modality

# Setup and train TOTALVI model
sca.models.TOTALVI.setup_anndata(
    adata_RNA,
    batch_key="batch",
    protein_expression_obsm_key="protein_expression"
)
arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
)
vae_ref = sca.models.TOTALVI(
    adata_RNA,
    **arches_params
)
vae_ref.train()

# Get latent representation
adata_RNA.obsm["X_scArches"] = vae_ref.get_latent_representation()
latent = adata_RNA.obsm['X_scArches'].copy()
np.savetxt(RNA_path + "/scArches.csv", latent, delimiter=',')

# 处理完成后保存 metadata.csv
adata_RNA.obs[['cell_type']].to_csv(RNA_path + "/metadata.csv")