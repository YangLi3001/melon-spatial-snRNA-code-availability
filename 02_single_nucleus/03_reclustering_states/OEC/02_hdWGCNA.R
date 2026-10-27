library(Seurat)
library(hdWGCNA)
library(WGCNA)

set.seed(12345)
enableWGCNAThreads(nThreads = 8)
data.oec <- readRDS("oec.rds")
DefaultAssay(data.oec) <- "SCT"

data.oec <- SetupForWGCNA(data.oec, wgcna_name = "OEC",
                          gene_select = "fraction", fraction = 0.05)
data.oec <- MetacellsByGroups(data.oec, group.by = "seurat_clusters",
                              k = 20, max_shared = 50, reduction = "umap",
                              ident.group = "seurat_clusters")
data.oec <- NormalizeMetacells(data.oec)
data.oec <- SetDatExpr(data.oec, assay = "SCT", layer = "data")
data.oec <- TestSoftPowers(data.oec, networkType = "signed")
data.oec <- ConstructNetwork(data.oec, soft_power = 4, tom_name = "OEC",
                             setDatExpr = FALSE, overwrite_tom = TRUE)

modules <- GetModules(data.oec)
write.csv(modules, "OEC_modules.csv", row.names = FALSE)
saveRDS(data.oec, "OEC_hdWGCNA.rds")
