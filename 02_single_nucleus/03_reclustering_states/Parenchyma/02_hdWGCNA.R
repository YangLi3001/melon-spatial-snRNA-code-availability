library(Seurat)
library(hdWGCNA)
library(WGCNA)

set.seed(12345)
enableWGCNAThreads(nThreads = 8)
data.BB <- readRDS("parenchyma.rds")
DefaultAssay(data.BB) <- "SCT"

data.BB <- SetupForWGCNA(data.BB, wgcna_name = "Pa",
                         gene_select = "fraction", fraction = 0.05)
data.BB <- MetacellsByGroups(data.BB, group.by = "seurat_clusters",
                             k = 20, max_shared = 30, reduction = "umap",
                             ident.group = "seurat_clusters")
data.BB <- NormalizeMetacells(data.BB)
data.BB <- SetDatExpr(data.BB, assay = "SCT", layer = "data")
data.BB <- TestSoftPowers(data.BB, networkType = "signed")
data.BB <- ConstructNetwork(data.BB, soft_power = 5, tom_name = "Pa",
                            setDatExpr = FALSE, overwrite_tom = TRUE)

modules <- GetModules(data.BB)
write.csv(modules, "parenchyma_modules.csv", row.names = FALSE)
saveRDS(data.BB, "parenchyma_hdWGCNA.rds")
