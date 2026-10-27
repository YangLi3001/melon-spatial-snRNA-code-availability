# GO enrichment with the melon gene-to-GO annotation.
# Usage: Rscript 02_GO_enrichment.R gene_list.txt go_annotation.tsv go_terms.tsv

suppressPackageStartupMessages(library(clusterProfiler))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Provide: gene_list, gene_to_GO table and GO term table")
}

genes <- unique(read.delim(args[[1]], header = FALSE,
                           stringsAsFactors = FALSE)[, 1])
gene_to_go <- read.delim(args[[2]], header = FALSE,
                         stringsAsFactors = FALSE)
colnames(gene_to_go)[1:2] <- c("gene_id", "ID")
go_terms <- read.delim(args[[3]], header = FALSE,
                       stringsAsFactors = FALSE)
colnames(go_terms)[1:3] <- c("ID", "Ontology", "Description")

go <- enricher(
  gene = genes,
  TERM2GENE = gene_to_go[c("ID", "gene_id")],
  TERM2NAME = go_terms[c("ID", "Description")],
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  maxGSSize = 500
)

write.table(as.data.frame(go), paste0(args[[1]], ".go.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
