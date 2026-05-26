#!/usr/bin/env Rscript
# edger.R — edgeR differential expression and TMM normalisation
#
# Usage:  Rscript edger.R <params.json>
#
# The params JSON must contain:
#   mode              : "unpaired" | "tmm"
#   counts_parquet    : path to count matrix (Parquet)
#   output_tsv        : output path
#   output_parquet    : output path
#
# For mode == "unpaired":
#   condition_a             : treatment condition label
#   condition_b             : reference condition label
#   columns_a               : list of sample column names for condition_a
#   columns_b               : list of sample column names for condition_b
#   comparison_name         : human-readable label (for logging)
#   library_sizes           : list of floats, or null (use column sums)
#   manual_dispersion_value : float, fallback dispersion for n=1 replicates
#   column_map              : named object R-col -> pipeline-col
#
# For mode == "tmm":
#   sample_columns    : ordered list of sample column names to include
#   library_sizes     : list of floats, or null (use column sums)
#   log               : logical, return log2-CPM (default TRUE)
#   prior_count       : prior count added before log (default 2)

suppressPackageStartupMessages({
  library(edgeR)
  library(arrow)
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# 1.  Read params
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript edger.R <params.json>", call. = FALSE)
}

p <- fromJSON(args[1], simplifyVector = FALSE)

mode           <- p$mode
output_tsv     <- p$output_tsv
output_parquet <- p$output_parquet

# ---------------------------------------------------------------------------
# 2.  Load counts
# ---------------------------------------------------------------------------
message("Reading counts from ", p$counts_parquet)
counts_df <- as.data.frame(read_parquet(p$counts_parquet))
rownames(counts_df) <- counts_df[[1]]
counts_df[[1]] <- NULL

# ---------------------------------------------------------------------------
# 3.  Dispatch
# ---------------------------------------------------------------------------
if (mode == "unpaired") {

  columns_a <- unlist(p$columns_a)
  columns_b <- unlist(p$columns_b)
  col_map   <- p$column_map
  n_a <- length(columns_a)
  n_b <- length(columns_b)

  all_cols  <- c(columns_a, columns_b)
  count_mat <- round(as.matrix(counts_df[, all_cols]))

  # Rename to X_0, X_1, … to avoid edgeR complaints about special chars
  renamed <- paste0("X_", seq_len(ncol(count_mat)) - 1)
  colnames(count_mat) <- renamed

  lib_sizes <- if (!is.null(p$library_sizes)) unlist(p$library_sizes) else colSums(count_mat)

  samples_df <- data.frame(
    group    = factor(c(rep("z", n_a), rep("x", n_b))),  # z=treatment, x=ref
    lib.size = lib_sizes,
    row.names = renamed
  )

  dge <- DGEList(counts = count_mat, samples = samples_df)
  dge <- calcNormFactors(dge)

  if (n_a == 1 && n_b == 1) {
    message(
      "Single replicate per condition — using manual dispersion ",
      p$manual_dispersion_value,
      ".\nConsider adding more replicates or using logFC only."
    )
    exact_tested <- exactTest(dge, dispersion = (p$manual_dispersion_value)^2)
  } else {
    disp         <- estimateDisp(dge, robust = TRUE)
    exact_tested <- exactTest(disp)
  }

  res    <- topTags(exact_tested, n = nrow(count_mat), sort.by = "none")
  df_res <- as.data.frame(res$table)
  df_res <- df_res[rownames(counts_df), ]     # restore original gene order

  for (r_col in names(col_map)) {
    pipeline_col <- col_map[[r_col]]
    if (r_col %in% colnames(df_res)) {
      colnames(df_res)[colnames(df_res) == r_col] <- pipeline_col
    }
  }

  dir.create(dirname(output_tsv),     recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(output_parquet), recursive = TRUE, showWarnings = FALSE)

  message("Writing TSV     -> ", output_tsv)
  write.table(df_res, file = output_tsv, sep = "\t", quote = FALSE, col.names = NA)
  message("Writing Parquet -> ", output_parquet)
  write_parquet(df_res, output_parquet)
  message("edgeR unpaired done.")

} else if (mode == "tmm") {

  sample_columns <- unlist(p$sample_columns)
  do_log         <- if (!is.null(p$log))         as.logical(p$log)         else TRUE
  prior_count    <- if (!is.null(p$prior_count)) as.numeric(p$prior_count) else 2

  count_mat <- round(as.matrix(counts_df[, sample_columns]))
  colnames(count_mat) <- sample_columns   # keep original names in output

  lib_sizes <- if (!is.null(p$library_sizes)) unlist(p$library_sizes) else colSums(count_mat)

  dge <- DGEList(counts = count_mat, lib.size = lib_sizes)
  dge <- calcNormFactors(dge, method = "TMM")

  message(
    "Applying TMM normalisation (log2=", do_log,
    ", prior_count=", prior_count, ")"
  )
  mat <- cpm(dge, log = do_log, prior.count = prior_count)

  # Output: genes x samples, gene_id as first column
  df_res <- as.data.frame(mat)
  df_res <- cbind(gene_id = rownames(df_res), df_res)
  rownames(df_res) <- NULL

  dir.create(dirname(output_tsv),     recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(output_parquet), recursive = TRUE, showWarnings = FALSE)

  message("Writing TSV     -> ", output_tsv)
  write.table(df_res, file = output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Writing Parquet -> ", output_parquet)
  write_parquet(df_res, output_parquet)
  message("edgeR TMM done.")

} else {
  stop(paste("Unknown mode:", mode), call. = FALSE)
}