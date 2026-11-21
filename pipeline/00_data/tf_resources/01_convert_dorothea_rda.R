args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  input_rda <- args[1]
  output_tsv <- args[2]
} else {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  script_path <- normalizePath(sub("--file=", "", cmd_args[grep("^--file=", cmd_args)]), winslash = "/", mustWork = TRUE)
  repo <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
  data_root <- Sys.getenv("CRC_DATA_ROOT", file.path(repo, "datasets"))
  input_rda <- file.path(data_root, "TF_resources", "data", "dorothea_hs_pancancer.rda")
  output_tsv <- file.path(data_root, "TF_resources", "data", "dorothea_hs_pancancer.tsv")
}

if (!file.exists(input_rda)) {
  stop(paste("Input RDA file not found:", input_rda))
}

object_names <- load(input_rda)
if (length(object_names) < 1) {
  stop("The RDA file holds no object")
}

obj <- get(object_names[1])
df <- as.data.frame(obj, stringsAsFactors = FALSE)

dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)
write.table(df, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Converted:", input_rda, "->", output_tsv, "\n")
