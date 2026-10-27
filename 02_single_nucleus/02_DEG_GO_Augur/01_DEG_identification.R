# Differential expression within cell types: P versus MO and C versus MO.
# Usage: Rscript 01_DEG_identification.R snRNA.rds [output_dir]
suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
rds_path <- if (length(args) >= 1) args[[1]] else "snRNA.rds"
output_dir <- if (length(args) >= 2) args[[2]] else "results/snRNA_DEG"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data.all <- readRDS(rds_path)
DefaultAssay(data.all) <- "RNA"

cluster_to_celltype <- c(
  "1" = "VT", "30" = "VT", "16" = "VT", "25" = "VT", "27" = "VT",
  "21" = "VT", "24" = "VT", "28" = "VT", "29" = "VT",
  "17" = "RT", "19" = "RT",
  "22" = "Ov", "23" = "Ov", "18" = "Ov", "4" = "Ov",
  "13" = "Ep", "14" = "Ep",
  "5" = "PC", "15" = "PC", "8" = "PC", "26" = "PC",
  "3" = "Pa", "11" = "Pa", "2" = "Pa", "12" = "Pa", "0" = "Pa",
  "6" = "Pa", "20" = "Pa", "10" = "Pa", "9" = "Pa", "7" = "Pa"
)
data.all$cell_type <- unname(cluster_to_celltype[as.character(data.all$seurat_clusters)])

comparisons <- c("P", "C")
control <- "MO"
min_cells <- 10
min_pct <- 0.1
padj_cutoff <- 0.05
logfc_cutoff <- 0.25

for (cell_type in c("VT", "RT", "Ov", "Ep", "PC", "Pa")) {
  union_genes <- character()
  for (sample_id in comparisons) {
    cells <- colnames(data.all)[
      !is.na(data.all$cell_type) & data.all$cell_type == cell_type &
        data.all$orig.ident %in% c(control, sample_id)
    ]
    if (length(cells) == 0) next
    obj <- subset(data.all, cells = cells)
    group_counts <- table(factor(obj$orig.ident, levels = c(control, sample_id)))
    if (any(group_counts < min_cells)) next
    Idents(obj) <- "orig.ident"
    markers <- FindMarkers(
      obj, ident.1 = sample_id, ident.2 = control,
      test.use = "wilcox", only.pos = FALSE,
      logfc.threshold = 0, min.pct = min_pct
    )
    if (nrow(markers) == 0) next
    markers$gene <- rownames(markers)
    logfc_cols <- intersect(c("avg_log2FC", "avg_logFC"), names(markers))
    stopifnot(length(logfc_cols) > 0)
    logfc_col <- logfc_cols[which.max(vapply(markers[logfc_cols], function(x) sum(!is.na(x)), numeric(1)))]
    deg <- markers[!is.na(markers$p_val_adj) & markers$p_val_adj < padj_cutoff &
                     !is.na(markers[[logfc_col]]) & abs(markers[[logfc_col]]) >= logfc_cutoff, ]
    prefix <- file.path(output_dir, paste(cell_type, sample_id, "vs", control, sep = "_"))
    write.csv(markers, paste0(prefix, "_all_markers.csv"), row.names = FALSE)
    write.csv(deg, paste0(prefix, "_DEGs.csv"), row.names = FALSE)
    union_genes <- union(union_genes, deg$gene)
  }
  writeLines(sort(union_genes), file.path(output_dir, paste0(cell_type, "_DEG_union.txt")))
}
