library(Seurat)
library(monocle)
library(Matrix)

# Cell-state annotations are added by 01_reclustering.R.
data.recluster <- readRDS("parenchyma.rds")
data.recluster <- subset(data.recluster, subset = orig.ident == "MO")
counts <- as(GetAssayData(data.recluster, assay = "SCT", slot = "counts"), "sparseMatrix")
pd <- new("AnnotatedDataFrame", data = data.recluster@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(counts),
                                                row.names = rownames(counts)))
cds <- newCellDataSet(counts, phenoData = pd, featureData = fd,
                      lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

diff <- differentialGeneTest(cds, fullModelFormulaStr = "~seurat_clusters", cores = 30)
deg <- subset(diff, qval < 0.01)
ordering_genes <- head(deg$gene_short_name[order(deg$qval)], 600)
cds <- setOrderingFilter(cds, ordering_genes)
HSMM <- reduceDimension(cds, max_components = 2, method = "DDRTree")
HSMM <- orderCells(HSMM, num_paths = 2)
HSMM <- orderCells(HSMM, root_state = 2)

BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 30)
BEAM_res <- BEAM_res[order(BEAM_res$qval), c("gene_short_name", "pval", "qval")]
write.csv(BEAM_res, "MO_BEAM_results.csv")
saveRDS(HSMM, "MO_monocle.rds")
pdf("MO_trajectory_heatmap.pdf")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_genes_branched_heatmap(HSMM[rownames(subset(BEAM_res, pval < 0.05)), ], branch_point = 1, num_clusters = 4, cores = 30, use_gene_short_name = TRUE, show_rownames = FALSE, return_heatmap = FALSE)
dev.off()
