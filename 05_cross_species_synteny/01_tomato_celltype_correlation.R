#!/usr/bin/env Rscript
# Usage: Rscript 01_tomato_celltype_correlation.R snRNA.rds tomato_markers.tsv homologs.tsv [output_dir]
# Tomato profiles are marker-derived pseudo-averages: 2^avg_log2FC for positive
# tissue markers and a baseline of 1 elsewhere.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(Rcpp)
})
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Provide the melon Seurat RDS, tomato markers and BLAST top-hit table")
RDS_FILE <- args[1]
TOMATO_MARKER_FILE <- args[2]
BLAST_TOP1_FILE <- args[3]
out_dir <- if (length(args) >= 4L) args[4] else "results/tomato_correlation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path, ...) {
  read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
             check.names = FALSE, ...)
}

normalize_tomato_id <- function(x) {
  x <- trimws(as.character(x))
  dot_count <- lengths(regmatches(x, gregexpr("\\.", x)))
  out <- x
  out[dot_count >= 2L] <- sub("\\.[0-9]+$", "", x[dot_count >= 2L])
  out
}

read_tomato_markers <- function() {
  x <- read_tsv(TOMATO_MARKER_FILE)
  names(x) <- make.names(names(x), unique = TRUE)
  x$gene <- normalize_tomato_id(x$gene_id)
  x$cluster <- tolower(trimws(as.character(x$cluster)))
  for (nm in intersect(c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), names(x))) {
    x[[nm]] <- as.numeric(x[[nm]])
  }
  x
}

filter_positive_markers <- function(x, min_fc = 0.25, min_pct_delta = 0.05,
                                    max_padj = 0.05) {
  keep <- is.finite(x$avg_log2FC) & x$avg_log2FC >= min_fc
  if ("p_val_adj" %in% names(x)) keep <- keep & x$p_val_adj <= max_padj
  if (all(c("pct.1", "pct.2") %in% names(x))) {
    keep <- keep & (x$pct.1 - x$pct.2) >= min_pct_delta
  }
  x[keep, , drop = FALSE]
}

read_blast <- function(path, max_hits = 1L, evalue_max = 1e-5,
                       pident_min = 25, align_min = 40) {
  cols <- c("melon_gene", "tomato_protein", "pident", "align_length", "mismatch",
            "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  x <- read.delim(path, header = FALSE, sep = "\t", quote = "", comment.char = "",
                  col.names = cols, stringsAsFactors = FALSE)
  x$melon_gene <- as.character(x$melon_gene)
  x$tomato_gene <- normalize_tomato_id(x$tomato_protein)
  x <- x[is.finite(x$evalue) & x$evalue <= evalue_max &
           is.finite(x$pident) & x$pident >= pident_min &
           is.finite(x$align_length) & x$align_length >= align_min, , drop = FALSE]
  x <- x[order(x$melon_gene, -x$bitscore, x$evalue, -x$pident), , drop = FALSE]
  key <- paste(x$melon_gene, x$tomato_gene, sep = "\r")
  x <- x[!duplicated(key), , drop = FALSE]
  x$hit_rank <- ave(x$bitscore, x$melon_gene, FUN = function(z) seq_along(z))
  x <- x[x$hit_rank <= max_hits, , drop = FALSE]
  max_bit <- ave(x$bitscore, x$melon_gene, FUN = max)
  x$map_weight <- x$bitscore / max_bit
  target_frequency <- table(x$tomato_gene)
  x$unique_target <- as.integer(target_frequency[x$tomato_gene]) == 1L
  rownames(x) <- NULL
  x
}

get_besthit_map <- function() read_blast(BLAST_TOP1_FILE, max_hits = 1L)

write_table <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

nperm <- as.integer(Sys.getenv("CORTABLE_NPERM", "1000"))
if (!is.finite(nperm) || nperm < 10L) stop("CORTABLE_NPERM must be at least 10")

celltype_by_cluster <- c(
  "0" = "Pa_1", "1" = "VT", "2" = "Pa_1", "3" = "Pa_2",
  "4" = "Ov", "5" = "PC", "6" = "Pa_1", "7" = "Pa_3",
  "8" = "PC", "9" = "Pa_3", "10" = "Pa_1", "11" = "Pa_2",
  "12" = "Pa_1", "13" = "Ep", "14" = "Ep", "15" = "PC",
  "16" = "VT", "17" = "RT", "18" = "Ov", "19" = "RT",
  "20" = "Pa_1", "21" = "VT", "22" = "Ov", "23" = "Ov",
  "24" = "VT", "25" = "VT", "26" = "PC", "27" = "VT",
  "28" = "VT", "29" = "VT", "30" = "VT"
)
celltype_levels <- c("Pa_1", "Pa_2", "Pa_3", "VT", "Ov", "PC", "Ep", "RT")

normalize_across_groups <- function(x) {
  x <- as.matrix(x)
  out <- sweep(x, 1, rowMeans(x), "/")
  out[!is.finite(out)] <- 0
  out
}

Rcpp::cppFunction('
Rcpp::NumericMatrix shuffle_rows_cpp(Rcpp::NumericMatrix x) {
  Rcpp::RNGScope scope;
  Rcpp::NumericMatrix out = Rcpp::clone(x);
  int nr = out.nrow();
  int nc = out.ncol();
  for (int i = 0; i < nr; ++i) {
    for (int j = nc - 1; j > 0; --j) {
      int k = static_cast<int>(R::unif_rand() * (j + 1));
      double tmp = out(i, j);
      out(i, j) = out(i, k);
      out(i, k) = tmp;
    }
  }
  return out;
}
')

cortable_cross <- function(ref, que, nperm = 1000L, seed = 1L) {
  common <- intersect(rownames(ref), rownames(que))
  ref <- as.matrix(ref[common, , drop = FALSE])
  que <- as.matrix(que[common, , drop = FALSE])
  observed <- suppressWarnings(stats::cor(que, ref, method = "spearman"))
  null_sum <- matrix(0, nrow(observed), ncol(observed), dimnames = dimnames(observed))
  exceed <- matrix(0L, nrow(observed), ncol(observed), dimnames = dimnames(observed))

  set.seed(seed)
  for (i in seq_len(nperm)) {
    z <- suppressWarnings(stats::cor(shuffle_rows_cpp(que), shuffle_rows_cpp(ref),
                                     method = "spearman"))
    z[!is.finite(z)] <- 0
    null_sum <- null_sum + z
    exceed <- exceed + (abs(z) >= abs(observed))
    if (i %% 100L == 0L || i == nperm) cat("permutations", i, "of", nperm, "\n")
  }

  p_original <- exceed / nperm
  p_plus1 <- (exceed + 1) / (nperm + 1)
  null_mean <- null_sum / nperm
  tab <- data.frame(
    melon_celltype = rep(rownames(observed), times = ncol(observed)),
    tomato_celltype = rep(colnames(observed), each = nrow(observed)),
    spearman = as.vector(observed),
    null_mean = as.vector(null_mean),
    spearman_minus_null = as.vector(observed - null_mean),
    p_permutation_original = as.vector(p_original),
    p_permutation_plus1 = as.vector(p_plus1),
    n_genes = length(common),
    stringsAsFactors = FALSE
  )
  tab$fdr_within_melon_celltype <- ave(
    tab$p_permutation_plus1, tab$melon_celltype,
    FUN = function(z) stats::p.adjust(z, method = "BH")
  )
  tab$fdr_global <- stats::p.adjust(tab$p_permutation_plus1, method = "BH")
  tab$tomato_rank <- ave(-tab$spearman, tab$melon_celltype,
                         FUN = function(z) rank(z, ties.method = "average"))
  tab$original_significance <- ifelse(
    tab$p_permutation_original < 0.01, "**",
    ifelse(tab$p_permutation_original < 0.05, "*", "")
  )
  tab$fdr_significance <- ifelse(
    tab$fdr_within_melon_celltype < 0.01, "**",
    ifelse(tab$fdr_within_melon_celltype < 0.05, "*", "")
  )
  list(correlation = observed, p_original = p_original,
       p_plus1 = p_plus1, null_mean = null_mean, table = tab)
}

cat("Reading melon Seurat object\n")
obj <- readRDS(RDS_FILE)
meta <- obj@meta.data
clusters <- if ("seurat_clusters" %in% names(meta)) {
  as.character(meta$seurat_clusters)
} else {
  as.character(SeuratObject::Idents(obj))
}
if (any(!clusters %in% names(celltype_by_cluster))) stop("Unmapped melon cluster detected")
melon_celltypes <- unname(celltype_by_cluster[clusters])
assay <- if ("SCT" %in% names(obj@assays)) "SCT" else SeuratObject::DefaultAssay(obj)

tomato <- filter_positive_markers(read_tomato_markers())
tomato_tissues <- sort(unique(tomato$cluster))
map_all <- get_besthit_map()
map <- map_all[map_all$unique_target, , drop = FALSE]
map <- map[map$melon_gene %in% rownames(obj[[assay]]) &
             map$tomato_gene %in% unique(tomato$gene), , drop = FALSE]
map <- map[!duplicated(map$melon_gene) & !duplicated(map$tomato_gene), , drop = FALSE]
map <- map[order(map$tomato_gene), , drop = FALSE]

obj$cortable_celltype <- factor(melon_celltypes, levels = celltype_levels)
cat("Calculating AverageExpression for", length(celltype_levels), "melon cell types and",
    nrow(map), "mapped genes\n")
avg_list <- Seurat::AverageExpression(
  obj, assays = assay, features = map$melon_gene,
  group.by = "cortable_celltype", slot = "data", verbose = FALSE
)
melon_avg <- as.matrix(avg_list[[assay]])
colnames(melon_avg) <- sub("^Pa-([123])$", "Pa_\\1", colnames(melon_avg))
missing_celltypes <- setdiff(celltype_levels, colnames(melon_avg))
if (length(missing_celltypes)) stop("Missing cell types: ", paste(missing_celltypes, collapse = ","))
melon_avg <- melon_avg[map$melon_gene, celltype_levels, drop = FALSE]
rownames(melon_avg) <- map$tomato_gene
melon_relative <- normalize_across_groups(melon_avg)

genes <- rownames(melon_relative)
tomato_fold <- matrix(1, nrow = length(genes), ncol = length(tomato_tissues),
                      dimnames = list(genes, tomato_tissues))
tm <- tomato[tomato$gene %in% genes & tomato$cluster %in% tomato_tissues,
             c("gene", "cluster", "avg_log2FC"), drop = FALSE]
tm <- aggregate(avg_log2FC ~ gene + cluster, tm, max)
tm$pseudo_average_fold <- 2^pmax(tm$avg_log2FC, 0)
tomato_fold[cbind(match(tm$gene, genes), match(tm$cluster, tomato_tissues))] <-
  tm$pseudo_average_fold
tomato_relative <- normalize_across_groups(tomato_fold)

write_table(data.frame(melon_gene = map$melon_gene, tomato_gene = map$tomato_gene),
            file.path(out_dir, "homolog_mapping.tsv"))
write_table(data.frame(gene = genes, melon_relative, check.names = FALSE),
            file.path(out_dir, "melon_relative_expression.tsv"))
write_table(data.frame(gene = genes, tomato_relative, check.names = FALSE),
            file.path(out_dir, "tomato_marker_derived_relative_expression.tsv"))

ans <- cortable_cross(tomato_relative, melon_relative, nperm = nperm,
                      seed = 2026072808L)
write_table(ans$table, file.path(out_dir, "celltype_correlation_results.tsv"))
write_table(data.frame(melon_celltype = rownames(ans$correlation), ans$correlation,
                       check.names = FALSE),
            file.path(out_dir, "celltype_spearman_matrix.tsv"))
