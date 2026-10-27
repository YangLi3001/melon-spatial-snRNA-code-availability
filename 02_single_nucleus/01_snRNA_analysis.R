# Seurat preprocessing and clustering.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

matrix_dirs <- c(
  MO = "data/MO",
  C  = "data/C",
  P  = "data/P",
  U  = "data/U"
)

read_filtered_matrix <- function(matrix_dir, sample_id) {
  counts <- readMM(file.path(matrix_dir, "matrix.mtx.gz"))
  features <- read.delim(file.path(matrix_dir, "features.tsv.gz"),
                         header = FALSE, stringsAsFactors = FALSE)[, 1]
  barcodes <- read.delim(file.path(matrix_dir, "barcodes.tsv.gz"),
                         header = FALSE, stringsAsFactors = FALSE)[, 1]
  rownames(counts) <- features
  colnames(counts) <- barcodes
  CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
}

snRNA <- Map(read_filtered_matrix, matrix_dirs, names(matrix_dirs))

snRNA_sct <- lapply(snRNA, function(x) {
  SCTransform(x, verbose = FALSE)
})

data.all <- merge(snRNA_sct[[1]], snRNA_sct[-1],
                  add.cell.ids = names(snRNA_sct))
DefaultAssay(data.all) <- "SCT"
VariableFeatures(data.all) <- unique(unlist(lapply(snRNA_sct, VariableFeatures)))

pc_dims <- 1:42
data.all <- RunPCA(data.all, verbose = FALSE)
data.all <- FindNeighbors(data.all, dims = pc_dims)

# Clustering parameter used in this workflow.
cluster_resolution <- 1.0
data.all <- FindClusters(data.all, resolution = cluster_resolution)

saveRDS(data.all, file = "snRNA.rds")
