library(Seurat)
library(Matrix)

# 1. File directory
rna_file <- "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSM4972212/GSM4972212_Citeseq_Human.GBM.R2_5_ND8.filtered.RNA.feature.bc.matrix.csv.gz"
adt_file <- "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSM4972212/GSM4972212_Citeseq_Human.GBM.R2_5_ND8.filtered.ADT.feature.bc.matrix.csv.gz"

# 2. Read csv
#    Assume the first column are feature (gene/ADT) name，the first row is barcode
M_rna_df <- read.csv(gzfile(rna_file), header = TRUE, row.names = 1, check.names = FALSE)
M_adt_df <- read.csv(gzfile(adt_file), header = TRUE, row.names = 1, check.names = FALSE)

# 3. Convert to a sparse matrix
M_rna <- as(Matrix(as.matrix(M_rna_df), sparse = TRUE), "dgCMatrix")
M_adt <- as(Matrix(as.matrix(M_adt_df), sparse = TRUE), "dgCMatrix")

# 4. Create SeuratObject
seu <- CreateSeuratObject(
  counts = M_rna,
  assay = "RNA",
  project = "GBM_CITEseq",
  min.cells = 3,
  min.features = 200
)

# 5. Add ADT as another assay
adt_assay <- CreateAssayObject(counts = M_adt)
seu[["ADT"]] <- adt_assay

# 6. Continue downstream analysis, such as standardization and integration
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

# 7. Multimodal analysis example
dim(Embeddings(seu, "pca.rna"))
# [1] 21589    50
dim(Embeddings(seu, "pca.adt")) 
# [1] 21589    37
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
path = "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSM4972212"
writeMM(seu@graphs$wsnn,paste0(path, "/Seurat_connectivities.mtx"))
writeMM(seu@graphs$wknn,paste0(path, "/Seurat_distance.mtx")) 


umap_df <- Embeddings(seu, reduction = "wnn.umap") %>% as.data.frame()
umap_df$cell_id <- rownames(umap_df)

path <- "/Users/zihang/Documents/CMML/CMML_ICA2/multimodal/GSM4972212/"
write.csv(umap_df, file = paste0(path, "Seurat_wnn_umap.csv"), row.names = FALSE)




