args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
repo <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
data_root <- Sys.getenv("CRC_DATA_ROOT", file.path(repo, "datasets"))
root <- file.path(data_root, "ImmuneDeconv")
data_dir <- file.path(root, "data")
results_dir <- file.path(root, "results")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("MCPcounter", quietly = TRUE)) {
  stop("The MCPcounter package is not installed.")
}

expr_file <- file.path(data_root, "TCGA", "data", "expression", "TCGA-BRCA.star_tpm.tsv.gz")
anno_file <- file.path(data_root, "TCGA", "results", "TCGA_BRCA_expression_gene_annotation_matched.tsv")

expr <- read.delim(gzfile(expr_file), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
anno <- read.delim(anno_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
expr$ensembl_gene_id <- sub("\\..*$", "", expr$Ensembl_ID)
expr <- merge(anno[, c("ensembl_gene_id", "gene_name")], expr, by = "ensembl_gene_id")
expr <- expr[!is.na(expr$gene_name) & expr$gene_name != "", ]

mat <- as.matrix(expr[, !(colnames(expr) %in% c("ensembl_gene_id", "gene_name", "Ensembl_ID")), drop = FALSE])
mode(mat) <- "numeric"
rownames(mat) <- expr$gene_name
mat <- rowsum(mat, group = rownames(mat), reorder = FALSE) / as.vector(table(rownames(mat)))

genes_file <- file.path(data_dir, "MCPcounter_genes.txt")
if (!file.exists(genes_file)) {
  stop("MCP-counter signature file not found: ", genes_file)
}
genes_signature <- read.table(genes_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)
scores <- MCPcounter::MCPcounter.estimate(mat, featuresType = "HUGO_symbols", genes = genes_signature)
scores_df <- data.frame(cell_population = rownames(scores), scores, check.names = FALSE)
write.table(scores_df, file.path(results_dir, "TCGA_BRCA_MCPcounter_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("MCP-counter scores:", nrow(scores), "populations x", ncol(scores), "samples\n")
