library(Seurat)
library(SeuratDisk)

# Cell-state annotations and UMAP are added by 01_reclustering.R.
data.recluster <- readRDS("parenchyma.rds")
stopifnot("celltype" %in% names(data.recluster@meta.data))
Idents(data.recluster) <- "celltype"
SaveH5Seurat(data.recluster, filename = "parenchyma.h5seurat", overwrite = TRUE)
Convert("parenchyma.h5seurat", dest = "h5ad", overwrite = TRUE)
