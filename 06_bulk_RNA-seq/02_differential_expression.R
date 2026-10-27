suppressPackageStartupMessages(library(DESeq2))

# Input: counts.csv from the StringTie workflow.
# Usage: Rscript 02_differential_expression.R counts.csv results
args <- commandArgs(trailingOnly = TRUE)
count_file <- if (length(args) >= 1L) args[1] else "counts.csv"
out_dir <- if (length(args) >= 2L) args[2] else "results/bulk_DEG"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mutant <- c("KO1", "KO2")
wild_type <- c("WT1", "WT2", "WT3")
samples <- c(mutant, wild_type)
lfc_cutoff <- 1
padj_cutoff <- 0.05

count_matrix <- read.csv(count_file, row.names = 1, check.names = FALSE)
count_matrix <- round(as.matrix(count_matrix))
storage.mode(count_matrix) <- "integer"
stopifnot(all(samples %in% colnames(count_matrix)))
sample_info <- data.frame(
  sample = samples,
  condition = factor(c(rep("mutant", length(mutant)),
                       rep("wild_type", length(wild_type))),
                     levels = c("wild_type", "mutant")),
  row.names = samples
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix[, samples, drop = FALSE],
  colData = sample_info,
  design = ~ condition
)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(
  dds, contrast = c("condition", "mutant", "wild_type"), alpha = padj_cutoff
))
res$gene_id <- rownames(res)
res$direction <- "none"
res$direction[!is.na(res$padj) & res$padj < padj_cutoff &
                res$log2FoldChange >= lfc_cutoff] <- "up"
res$direction[!is.na(res$padj) & res$padj < padj_cutoff &
                res$log2FoldChange <= -lfc_cutoff] <- "down"
res <- res[order(res$padj, -res$log2FoldChange, na.last = TRUE), ]
res <- res[c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
             "pvalue", "padj", "direction")]
write.table(res, file.path(out_dir, "DESeq2_all_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res[res$direction != "none", ],
            file.path(out_dir, "DESeq2_significant_DEGs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
normalized <- data.frame(gene_id = rownames(dds), counts(dds, normalized = TRUE),
                         check.names = FALSE)
write.table(normalized, file.path(out_dir, "normalized_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_info, file.path(out_dir, "sample_information.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
