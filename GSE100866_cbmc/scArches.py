import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import scarches as sca


# RNA_path = "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/"
RNA_path = "/public/home/chenjiaminggroup/shenibdproject/multimodal/GSE100866_cbmc/"
adata_RNA = sc.read_h5ad(RNA_path + "cbmc_rna.h5ad")
adata_adt = sc.read_h5ad(RNA_path + "cbmc_adt.h5ad")

adata_RNA.obs.rename(columns={"rna_annotations": "cell_type"}, inplace=True)
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
adata_RNA.obsm["protein_expression"] = adata_adt.X.toarray()  
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



