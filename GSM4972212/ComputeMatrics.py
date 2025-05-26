import scib
import anndata as ad
import pandas as pd
import numpy as np
from multiprocessing import Pool 
from scipy.io import mmread
from scipy.sparse import csr_matrix
import muon
import os
import scarches as sca
import scanpy as sc
from scib_metrics.benchmark import Benchmarker
import scib_metrics
from typing import Optional
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

# dt_list = ['Dataset2','Dataset3','Dataset4','Dataset5','Dataset6','Dataset8','Dataset9',
#            'Dataset10','Dataset13','Dataset17','Dataset24','Dataset26']

dt_list = ['GSM4972212']
path = '/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/'

# for dataset in dt_list:
#     latent = pd.read_csv(path + dataset + '/scArches.csv', header = None)
#     latent.drop(0, axis=1)
#     latent.to_csv(path + dataset + '/scArches.csv', header = None)



def count_metrics(dataset):
    # dataset = 'GSM4972212'
    data_path = path + dataset
    RNA_path = data_path
    rna_data = pd.read_csv(RNA_path + '/GSM4972212_Citeseq_Human.GBM.R2_5_ND8.filtered.RNA.feature.bc.matrix.csv.gz', index_col=0).T
    cell_names = rna_data.index  # Barcodes
    gene_names = rna_data.columns  # Features
    M_rna = csr_matrix(rna_data.values)
    # Read annotation data (if needed for downstream analysis)
    annot = pd.read_csv(
        RNA_path + '/GSM4972212_annot.Citeseq_Human.GBM.R2_5_ND8.csv.gz',
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
    # Read ADT data
    adt_data = pd.read_csv(RNA_path + '/GSM4972212_Citeseq_Human.GBM.R2_5_ND8.filtered.ADT.feature.bc.matrix.csv.gz', index_col=0).T
    # adt_features = adt_data.columns  # ADT features
    # M_adt = csr_matrix(adt_data.values)
    adata_RNA.obsm["protein_expression"] = adt_data
    adata = adata_RNA.copy()
    # read ADT data
    metadata = pd.read_csv(data_path + "/metadata.csv")
    metadata['cell_type'].index = adata.obs_names
    adata.obs['cell_type'] = metadata['cell_type'].astype('category')
    if np.where(adata.obs["cell_type"].isna())[0].shape[0]!=0:
        adata.obs["cell_type"] = adata.obs["cell_type"].cat.add_categories(['NaN'])
        adata.obs["cell_type"][np.where(adata.obs["cell_type"].isna())[0]] = 'NaN'
    adata.obs['batch'] = ['batch1'] * adata.shape[0]
    result = pd.DataFrame()
    metrics_list = []
    # method_list = ['TotalVI','scArches','Multigrate','SCOIT','MOJITOO','DeepMaps','scVAEIT']
    # method_list = ['scArches']
    method_list = ['scArches',"PCA"]
    for method in method_list:
        if os.path.exists(path + dataset + '/' + method + '.csv'):
            # method = 'scArches'
            # method = 'PCA'
            latent = pd.read_csv(path + dataset + '/' + method + '.csv', header = None)
            latent.index = adata.obs_names
            adata.obsm[method] = latent
            latentPCA = pd.read_csv(path + dataset + '/' + method + '.csv', header = None)
            latentPCA.index = adata.obs_names
            adata.obsm[method] = latentPCA
            sc.pp.neighbors(adata, use_rep=method)
            sc.tl.umap(adata)
            sc.tl.leiden(adata, key_added=f"cluster_{method}")
            # 绘制并优化 UMAP 图
            sc.pl.umap(
                adata, 
                color='cell_type', 
                size=10, 
                alpha=0.7, 
                frameon=False, 
                title=f'UMAP for {method}', 
                legend_loc='right margin', 
                show=False  # 防止在 notebook 中显示图像
            )
            save_path = path + dataset + '/' + method + 'umap.png'
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close()
            scib.metrics.cluster_optimal_resolution(adata, cluster_key=f"cluster_{method}", label_key="cell_type")
            ari = scib.metrics.ari(adata, cluster_key=f"cluster_{method}", label_key="cell_type")
            iso_asw = scib.metrics.isolated_labels_asw(adata, label_key="cell_type", batch_key='batch', embed=method,  verbose = False)
            nmi = scib.metrics.nmi(adata, cluster_key=f"cluster_{method}", label_key="cell_type")
            clisi = scib.metrics.clisi_graph(adata, label_key="cell_type",use_rep=method, type_='embed')
            sht = scib.metrics.silhouette(adata, label_key="cell_type", embed=method, metric='euclidean', scale=True)
            metrics_list.append([ari, iso_asw, nmi, clisi, sht, method])
    # method_list = ['Seurat','CiteFuse']
    method_list = ['Seurat']
    for method in method_list:
        if os.path.exists(path+ dataset + '/'+ method + '_connectivities.mtx'):
            # method = 'Seurat'
            con = mmread(path+ dataset + '/'+ method + '_connectivities.mtx')
            dis = mmread(path+ dataset + '/'+ method + '_distance.mtx')
            adata.uns['neighbors'] = {'connectivities_key': 'connectivities', 'distances_key': 'distances', 
                                      'params': {'n_neighbors': 20, 'method': 'umap', 'random_state': 0, 
                                                 'metric': 'euclidean'}}
            adata.uns['neighbors']['distance'] = csr_matrix(dis)
            adata.uns['neighbors']['connectivities'] = csr_matrix(con)
            adata.obsp['distance'] = csr_matrix(dis)
            adata.obsp['connectivities'] = csr_matrix(con)
            sc.tl.umap(adata, n_components=20)
            scib.metrics.cluster_optimal_resolution(adata, cluster_key=f"cluster_{method}", label_key="cell_type")
            ari = scib.metrics.ari(adata, cluster_key=f"cluster_{method}", label_key="cell_type")
            iso_asw = scib.metrics.isolated_labels_asw(adata, label_key="cell_type", batch_key='batch', embed="X_umap",  verbose = False)
            nmi = scib.metrics.nmi(adata, cluster_key=f"cluster_{method}", label_key="cell_type")
            clisi = scib.metrics.clisi_graph(adata, label_key="cell_type", type_='knn')
            sht = scib.metrics.silhouette(adata, label_key="cell_type", embed="X_umap", metric='euclidean', scale=True)
            metrics_list.append([ari, iso_asw, nmi, clisi, sht, method])
            umap_df = pd.read_csv(path + dataset + '/Seurat_wnn_umap.csv')
            umap_df.set_index('cell_id', inplace=True)
            # Verify index alignment
            if not umap_df.index.equals(adata.obs.index):
                common_cells = umap_df.index.intersection(adata.obs.index)
                if len(common_cells) == 0:
                    raise ValueError(f"No common cell IDs for {method}!")
                umap_df = umap_df.loc[common_cells]
                annotations_df = annotations_df.loc[common_cells]
                adata_subset = adata[common_cells].copy()
                print(f"Using {len(common_cells)} common cells for {method}.")
            else:
                adata_subset = adata
            adata_subset.obsm['X_umap'] = umap_df[['wnnUMAP_1', 'wnnUMAP_2']].values
            # 绘制并优化 UMAP 图
            sc.pl.umap(
                adata_subset, 
                color='cell_type', 
                size=10, 
                alpha=0.7, 
                frameon=False, 
                title=f'UMAP for {method}', 
                legend_loc='right margin', 
                show=False  # 防止在 notebook 中显示图像
            )
            save_path = path + dataset + '/' + method + 'umap.png'
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close()
    df = pd.DataFrame(metrics_list,columns = ['ARI','Isolated_Labels_ASW','NMI','cLISI_Graph','Silhouette','method'])
    result = pd.concat([result,df])
    result['Dataset'] = dataset
    result.to_csv(path + dataset + "/mc_result.csv",index = False)
    print(dataset)



for dataset in dt_list:
    count_metrics(dataset)


# pool = Pool(10)
# results = pool.map(count_metrics, dt_list)
# pool.close()
# pool.join()


