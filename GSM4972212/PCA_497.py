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

# Optional: Read annotation data (if needed for downstream analysis)
annot = pd.read_csv(
    RNA_path + 'GSM4972212_annot.Citeseq_Human.GBM.R2_5_ND8.csv.gz', 
    sep='\t',  # 根据实际分隔符调整
    index_col=False
)
# metadata = pd.read_csv(
#     RNA_path + 'GSM4972212_annot.Citeseq_Human.GBM.R2_5_ND8.csv.gz', 
#     sep='\t',  # 根据实际分隔符调整
#     index_col=False
# )

# Create AnnData object with RNA data
adata_RNA = ad.AnnData(M_rna, obs=pd.DataFrame(index=cell_names), var=pd.DataFrame(index=gene_names))
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
# # 创建一个全 0 的 DataFrame，行数匹配 adata_RNA 中的细胞数量，列数匹配 ADT 特征数量
# adt_data_full = pd.DataFrame(
#     np.zeros((len(adata_RNA.obs_names), len(adt_features))),
#     index=adata_RNA.obs_names,
#     columns=adt_features
# )
# adt_cell_names = adt_data.index
# # 找到 M_adt 的细胞名称与 adata_RNA 的细胞名称的交集
# common_cells = adt_cell_names[adt_cell_names.isin(adata_RNA.obs_names)]
# # 获取 M_adt 中对应于共有细胞的行
# idx = adt_cell_names.isin(common_cells)
# M_adt_common = M_adt.toarray()[idx]
# # 将数据填充到 adt_data_full 的对应行
# adt_data_full.loc[common_cells] = M_adt_common
# # 如果需要，将结果赋值给 adata_RNA.obsm
# adata_RNA.obsm["protein_expression"] = adt_data_full

# Add batch information
modality = ['batch1'] * adata_RNA.shape[0]
adata_RNA.obs['batch'] = modality

# 合并 RNA 和 ADT 数据
rna_dense = adata_RNA.X.toarray()  # 将 RNA 数据转换为稠密矩阵
adt_dense = adata_RNA.obsm["protein_expression"].values  # 获取 ADT 数据
combined_data = np.hstack((rna_dense, adt_dense))  # 水平拼接 RNA 和 ADT 数据

# 创建一个新的 AnnData 对象用于 PCA
adata_combined = ad.AnnData(
    combined_data,
    obs=adata_RNA.obs,
    var=pd.DataFrame(index=adata_RNA.var.index.tolist() + adata_RNA.obsm["protein_expression"].columns.tolist())
)

# 进行 PCA 降维
sc.pp.pca(adata_combined, n_comps=50)  # 设置 PCA 维度，可根据需要调整

# 获取 PCA 的 latent representation
latent = adata_combined.obsm['X_pca'].copy()

# 保存为 CSV 文件
np.savetxt(RNA_path + "/PCA.csv", latent, delimiter=',')

