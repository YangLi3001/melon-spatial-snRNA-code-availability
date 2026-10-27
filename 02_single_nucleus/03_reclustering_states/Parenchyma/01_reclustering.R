library(Seurat)

# Input: the full single-nucleus Seurat object.
data.all <- readRDS("snRNA.rds")
data.all$old.cluster <- data.all$seurat_clusters
Idents(data.all) <- "seurat_clusters"
sce <- subset(data.all, idents = c(12, 2, 6, 20, 0, 10))
DefaultAssay(sce) <- "SCT"
sce <- RunPCA(sce, verbose = FALSE)
sce <- FindNeighbors(sce, dims = 1:40, k.param = 30)
sce <- FindClusters(sce, resolution = 0.4)
sce <- RunUMAP(sce, dims = 1:40, min.dist = 1)

sce <- PrepSCTFindMarkers(sce, assay = "SCT", verbose = TRUE)
markers <- FindAllMarkers(sce, only.pos = FALSE, min.pct = 0.25,
                          logfc.threshold = 0.1, test.use = "wilcox")
write.table(markers, "parenchyma_markers.tsv", sep = "\t", quote = FALSE)

# Vegetative, initial, morphogenetic and dividing states.
state_by_cluster <- c("0" = "Mor", "1" = "Ini", "2" = "Veg", "3" = "Div",
                      "4" = "Div", "5" = "Veg", "6" = "Div", "7" = "Mor",
                      "8" = "Div", "9" = "Mor")
stopifnot(all(as.character(sce$seurat_clusters) %in% names(state_by_cluster)))
sce$celltype <- factor(unname(state_by_cluster[as.character(sce$seurat_clusters)]),
                       levels = c("Veg", "Ini", "Mor", "Div"))
saveRDS(sce, "parenchyma.rds")
