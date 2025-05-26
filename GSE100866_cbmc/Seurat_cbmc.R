library(SeuratData)
library(Seurat)
library(SeuratDisk)
# InstallData(ds = 'cbmc')
# InstallData(ds = 'thp1.eccite') # GSE153056 Human ECCITE-seq

data("cbmc") # GSE100866 cord blood mononuclear cells
cbmc <- UpdateSeuratObject(cbmc)
cbmc
# saveRDS(cbmc, file = "cbmc_updated.rds")
SaveH5Seurat(cbmc, filename = "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/cbmc.h5Seurat", overwrite = TRUE)
Convert("/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/cbmc.h5Seurat", dest = "h5ad", overwrite = TRUE)

DefaultAssay(cbmc) <- "RNA"
cbmc_rna <- subset(cbmc, features = rownames(cbmc[["RNA"]]))
str(cbmc_rna)
SaveH5Seurat(cbmc_rna, filename = "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/cbmc_rna.h5Seurat", overwrite = TRUE)
Convert("/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/cbmc_rna.h5Seurat", dest = "h5ad", overwrite = TRUE)

DefaultAssay(cbmc) <- "ADT"
cbmc_adt <- subset(cbmc, features = rownames(cbmc[["ADT"]]))
str(cbmc_adt)
SaveH5Seurat(cbmc_adt, filename = "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/cbmc_adt.h5Seurat", overwrite = TRUE)
Convert("/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/cbmc_adt.h5Seurat", dest = "h5ad", overwrite = TRUE)


meta_data <- cbmc@meta.data
cell_type_data <- data.frame(
  cell_id = rownames(meta_data),  # 提取细胞条形码作为索引
  cell_type = meta_data$rna_annotation  # 提取 cell_type 列
)
write.csv(cell_type_data, 
          file = "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/cell_type.csv", 
          row.names = FALSE)  # 不保存行名，因为 cell_id 已经作为一列

seu = cbmc
# Continue downstream analysis, such as standardization and integration
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, reduction.name = "pca.rna")

DefaultAssay(seu) <- "ADT"
seu <- NormalizeData(seu, normalization.method = "CLR", margin = 2)
seu <- FindVariableFeatures(seu, nfeatures = length(adt.idx))
seu <- ScaleData(seu)
seu <- RunPCA(seu, reduction.name = "pca.adt")

# Multimodal analysis example
dim(Embeddings(seu, "pca.rna"))
# [1] 8617   50
dim(Embeddings(seu, "pca.adt")) 
# [1] 8617    9
seu <- FindMultiModalNeighbors(
  seu,
  reduction.list = list("pca.rna", "pca.adt"),
  dims.list = list(1:30, 1:9),
  modality.weight.name = "RNA_Weights"  # 可选，记录各模态权重
)

# 8. UMAP 
seu <- RunUMAP(seu,
               nn.name = "weighted.nn",
               reduction.name = "wnn.umap",
               reduction.key  = "wnnUMAP_"
)

# 9. Clustering based on WNN graphs
seu <- FindClusters(seu,
                    graph.name = "wsnn",
                    algorithm  = 3,
                    resolution = 0.5,
                    verbose    = FALSE
)

str(seu@graphs$wsnn)
str(seu@graphs$wknn)
path = "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/"
writeMM(seu@graphs$wsnn,paste0(path, "/Seurat_connectivities.mtx"))
writeMM(seu@graphs$wknn,paste0(path, "/Seurat_distance.mtx")) 


DimPlot(
  seu,
  reduction = "wnn.umap",
  group.by  = "rna_annotations",
  label     = TRUE,            
  repel     = TRUE             
) + 
  ggtitle("WNN-UMAP: Seurat Clusters") +
  theme_minimal()



umap_df <- Embeddings(seu, reduction = "wnn.umap") %>% as.data.frame()
umap_df$cell_id <- rownames(umap_df)

path <- "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE100866_cbmc/"
write.csv(umap_df, file = paste0(path, "Seurat_wnn_umap.csv"), row.names = FALSE)

