library(SeuratDisk)
library(Seurat)
# Load the h5Seurat file
data <- LoadH5Seurat("~/Downloads/cux2+.h5Seurat",assays = "counts", reductions = "umap", graphs = F, neighbors = F, images = F)
DimPlot(data, reduction = "umap")
