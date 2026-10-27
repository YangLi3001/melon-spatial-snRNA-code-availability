library(Seurat)

# Input: the full single-nucleus Seurat object.
data.all <- readRDS("snRNA.rds")
data.all$old.cluster <- data.all$seurat_clusters
Idents(data.all) <- "seurat_clusters"
sce <- subset(data.all, idents = c(3, 11))
DefaultAssay(sce) <- "SCT"
sce <- RunPCA(sce, verbose = FALSE)
sce <- FindNeighbors(sce, dims = 1:20, k.param = 20)
sce <- FindClusters(sce, resolution = 0.3)
sce <- RunUMAP(sce, dims = 1:20, min.dist = 1)

sce <- PrepSCTFindMarkers(sce, assay = "SCT", verbose = TRUE)
markers <- FindAllMarkers(sce, only.pos = FALSE, min.pct = 0.25,
                          logfc.threshold = 0.585, test.use = "wilcox")
write.table(markers, "oec_markers.tsv", sep = "\t", quote = FALSE)
saveRDS(sce, "oec.rds")
