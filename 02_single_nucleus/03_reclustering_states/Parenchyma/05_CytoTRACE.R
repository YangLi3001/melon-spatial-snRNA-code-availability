library(Seurat)
library(CytoTRACE)

subdata <- readRDS("parenchyma.rds")
subdata <- subset(subdata, subset = orig.ident %in% c("MO", "C", "P"))
expr <- as.matrix(GetAssayData(subdata, assay = "RNA", slot = "counts"))

results <- CytoTRACE(expr)
saveRDS(results, "parenchyma_CytoTRACE.rds")

phenotype <- as.character(subdata$seurat_clusters)
names(phenotype) <- colnames(subdata)
umap <- Embeddings(subdata, reduction = "umap")
plotCytoTRACE(results, phenotype = phenotype, emb = umap)
