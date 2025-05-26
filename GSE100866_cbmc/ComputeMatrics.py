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
# RNA_path = "/public/home/chenjiaminggroup/shenibdproject/multimodal/Sca1pos_CITEseq/"
dt_list = ["GSE100866_cbmc"]
path = '/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/'

# for dataset in dt_list:
#     latent = pd.read_csv(path + dataset + '/scArches.csv', header = None)
#     latent.drop(0, axis=1)
#     latent.to_csv(path + dataset + '/scArches.csv', header = None)



def count_metrics(dataset):
    # dataset = 'GSE100866_cbmc'
    RNA_path = path + dataset
    adata_RNA = sc.read_h5ad(RNA_path + "/cbmc_rna.h5ad")
    adata_adt = sc.read_h5ad(RNA_path + "/cbmc_adt.h5ad")
    cell_type_df = pd.read_csv(RNA_path + "/cell_type.csv")
    cell_type_df.set_index('cell_id', inplace=True)
    adata_RNA.obs['cell_type'] = cell_type_df['cell_type']
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
    adata = adata_RNA.copy()
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
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
            if "neighbors" in adata.uns:
                del adata.uns["neighbors"]
                adata.obsp.pop("connectivities", None)
                adata.obsp.pop("distances", None)
            sc.pp.neighbors(adata, use_rep=method)
            sc.tl.umap(adata)
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
            sc.tl.leiden(adata, key_added=f"cluster_{method}")
            # sc.pl.umap(adata, color=f"cluster_{method}", show=False)
            # plt.savefig(path + dataset + '/' + method + 'cluster_umap.png')
            # plt.close()
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
            sc.tl.umap(adata)
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
            sc.pl.umap(
                adata_subset, 
                color='cell_type', 
                size=10, 
                alpha=0.7, 
                frameon=False, 
                title=f'UMAP for {method}', 
                legend_loc='right margin', 
                show=False  
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

