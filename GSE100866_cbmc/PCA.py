import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad

RNA_path = "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/"
adata_RNA = sc.read_h5ad("/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/cbmc_rna.h5ad")
adata_adt = sc.read_h5ad("/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/cbmc_adt.h5ad")

adata_RNA.obs.rename(columns={"rna_annotations": "cell_type"}, inplace=True)
adata_RNA.obsm["protein_expression"] = adata_adt.X.toarray()  


# Add batch information
modality = ['batch1'] * adata_RNA.shape[0]
adata_RNA.obs['batch'] = modality

# 合并 RNA 和 ADT 数据
rna_dense = adata_RNA.X.toarray()  # 将 RNA 数据转换为稠密矩阵
adt_dense = adata_RNA.obsm["protein_expression"] # 获取 ADT 数据
combined_data = np.hstack((rna_dense, adt_dense))  # 水平拼接 RNA 和 ADT 数据

# 构造变量名（特征名）
var_names = adata_RNA.var_names.tolist() + adata_adt.var_names.tolist()
var_df = pd.DataFrame(index=var_names)

# 构建 AnnData 对象
adata_combined = ad.AnnData(
    X=combined_data,
    obs=adata_RNA.obs.copy(),
    var=var_df
)

# 进行 PCA 降维
sc.pp.pca(adata_combined, n_comps=50)  # 设置 PCA 维度，可根据需要调整

# 获取 PCA 的 latent representation
latent = adata_combined.obsm['X_pca'].copy()

# 保存为 CSV 文件
np.savetxt(RNA_path + "/PCA.csv", latent, delimiter=',')




