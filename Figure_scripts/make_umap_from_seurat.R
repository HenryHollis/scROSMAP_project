library(SeuratDisk)
library(Seurat)
#convert to h5ad files:
Convert("~/Downloads/cux2+.h5Seurat", dest = "h5ad", assay = "RNA")
Convert("~/Downloads/cux2-.h5Seurat", dest = "h5ad", assay = "RNA")
# Load the h5Seurat file
data <- LoadH5Seurat("~/Downloads/cux2+.h5Seurat",assays = "counts", reductions = "umap", graphs = F, neighbors = F, images = F)
DimPlot(data, reduction = "umap")
