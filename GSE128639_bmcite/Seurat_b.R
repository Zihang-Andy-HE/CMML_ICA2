library(SeuratData)
library(Seurat)
library(SeuratDisk)
# InstallData(ds = 'cbmc') # GSE100866 cord blood mononuclear cells
# InstallData(ds = 'thp1.eccite') # GSE153056 Human ECCITE-seq Human bone marrow mononuclear cells
# options(timeout = 600)
InstallData(ds = "bmcite") # GSE128639 Human bone marrow mononuclear cells

data("bmcite") 
bmcite <- UpdateSeuratObject(bmcite)
bmcite


# DefaultAssay(bmcite) <- "RNA"
# bmcite_rna <- subset(bmcite, features = rownames(bmcite[["RNA"]]))
# str(bmcite_rna)
# SaveH5Seurat(bmcite_rna, filename = "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE128639_bmcite/bmcite_rna.h5Seurat", overwrite = TRUE)
# Convert("/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE128639_bmcite/bmcite_rna.h5Seurat", dest = "h5ad", overwrite = TRUE)
# 
# DefaultAssay(bmcite) <- "ADT"
# bmcite_adt <- subset(bmcite, features = rownames(bmcite[["ADT"]]))
# bmcite_adt[["ADT"]]@data <- bmcite_adt[["ADT"]]@counts
# str(bmcite_adt)
# SaveH5Seurat(bmcite_adt, filename = "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE128639_bmcite/bmcite_adt.h5Seurat", overwrite = TRUE)
# Convert("/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE128639_bmcite/bmcite_adt.h5Seurat", dest = "h5ad", slot = "counts", overwrite = TRUE)


seu = bmcite
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
# [1] 30672    50
dim(Embeddings(seu, "pca.adt")) 
# [1] 30672    24
seu <- FindMultiModalNeighbors(
  seu,
  reduction.list = list("pca.rna", "pca.adt"),
  dims.list = list(1:30, 1:10),
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
path = "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE128639_bmcite/"
writeMM(seu@graphs$wsnn,paste0(path, "Seurat_connectivities.mtx"))
writeMM(seu@graphs$wknn,paste0(path, "Seurat_distance.mtx")) 




DimPlot(
  seu,
  reduction = "wnn.umap",
  group.by  = "celltype.l2",
  label     = TRUE,            
  repel     = TRUE            
) + 
  ggtitle("WNN-UMAP: Seurat Clusters") +
  theme_minimal()



umap_df <- Embeddings(seu, reduction = "wnn.umap") %>% as.data.frame()
umap_df$cell_id <- rownames(umap_df)

path <- "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSE128639_bmcite/"
write.csv(umap_df, file = paste0(path, "Seurat_wnn_umap.csv"), row.names = FALSE)

# Dataset, Cell, RNA, ADT, Description
# GSM4972212, 21589, 14131, 268, Profiling of glioblastoma-associated myeloid cells human
# GSE100866, 8617, 20501, 10, Cord blood mononuclear cells mouse
# GSE128639, 30672, 17009, 25, Human bone marrow mononuclear cells human

