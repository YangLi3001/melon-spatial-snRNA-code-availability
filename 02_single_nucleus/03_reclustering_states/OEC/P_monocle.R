library(Seurat)
library(monocle)
library(Matrix)

obj <- readRDS("oec.rds")
obj$second_cluster <- obj$seurat_clusters
Idents(obj) <- "seurat_clusters"
obj <- subset(obj, idents = c(0, 3, 5, 6))
state_by_cluster <- c("0" = "Biosynthetic", "3" = "Biosynthetic",
                      "5" = "Transport", "6" = "Primed")
obj$celltype <- unname(state_by_cluster[as.character(obj$seurat_clusters)])
obj <- subset(obj, subset = orig.ident == "P")
counts <- as(GetAssayData(obj, assay = "SCT", slot = "counts"), "sparseMatrix")
pd <- new("AnnotatedDataFrame", data = obj@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(counts),
                                                row.names = rownames(counts)))
cds <- newCellDataSet(counts, phenoData = pd, featureData = fd,
                      lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

diff <- differentialGeneTest(cds, fullModelFormulaStr = "~celltype", cores = 48)
deg <- subset(diff, qval < 0.01)
ordering_genes <- head(deg$gene_short_name[order(deg$qval)], 500)
cds <- setOrderingFilter(cds, ordering_genes)
HSMM <- reduceDimension(cds, max_components = 2, method = "DDRTree")
HSMM <- orderCells(HSMM)

pseudotime_genes <- differentialGeneTest(HSMM, fullModelFormulaStr = "~sm.ns(Pseudotime)")
significant_genes <- rownames(subset(pseudotime_genes, qval < 0.01))
write.csv(pseudotime_genes, "P_pseudotime_gene_test.csv")
saveRDS(HSMM, "P_OEC_monocle.rds")
pdf("P_OEC_trajectory_heatmap.pdf")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_pseudotime_heatmap(HSMM[significant_genes, ], num_clusters = 4, show_rownames = FALSE, return_heatmap = FALSE)
dev.off()
