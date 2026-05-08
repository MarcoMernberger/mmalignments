#!/usr/bin/env Rscript
# deseq2.R — DESeq2 differential expression analysis
#
# Usage:  Rscript deseq2.R <params.json>
#
# The params JSON must contain:
#   mode              : "unpaired" | "timeseries"
#   counts_parquet    : path to input count matrix (Parquet)
#   column_map        : named object mapping R column names -> pipeline column names
#   output_tsv        : path for TSV output
#   output_parquet    : path for Parquet output
#
# For mode == "unpaired":
#   condition_a       : treatment condition label
#   condition_b       : reference condition label
#   model_conditions  : named list: condition -> [sample columns]
#                       (may include extra conditions for variance estimation)
#
# For mode == "timeseries":
#   sample_columns    : ordered list of sample column names
#   factors           : named list: factor_name -> [per-sample values]
#   formula           : full model formula, e.g. "~ condition + time"
#   reduced           : reduced formula for LRT, e.g. "~ condition"

suppressPackageStartupMessages({
  library(DESeq2)
  library(arrow)      # Parquet I/O
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# 1.  Read params
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript deseq2.R <params.json>", call. = FALSE)
}

p <- fromJSON(args[1], simplifyVector = FALSE)

mode           <- p$mode
col_map        <- p$column_map          # named list R-col -> pipeline-col
output_tsv     <- p$output_tsv
output_parquet <- p$output_parquet

# ---------------------------------------------------------------------------
# 2.  Load count matrix
# ---------------------------------------------------------------------------
message("Reading counts from ", p$counts_parquet)
counts_df <- as.data.frame(read_parquet(p$counts_parquet))
rownames(counts_df) <- counts_df[[1]]           # first column is gene index
counts_df[[1]] <- NULL

# ---------------------------------------------------------------------------
# 3.  Build colData and select columns
# ---------------------------------------------------------------------------
if (mode == "unpaired") {

  model_conditions <- p$model_conditions          # named list cond -> [cols]
  condition_a <- p$condition_a
  condition_b <- p$condition_b

  # Build sample metadata frame
  all_samples <- c()
  all_labels  <- c()
  for (cond in names(model_conditions)) {
    cols <- unlist(model_conditions[[cond]])
    all_samples <- c(all_samples, cols)
    all_labels  <- c(all_labels, rep(cond, length(cols)))
  }

  col_data <- data.frame(
    condition = factor(all_labels, levels = c(condition_a, condition_b,
                       setdiff(names(model_conditions), c(condition_a, condition_b)))),
    row.names = all_samples
  )
  count_mat <- round(as.matrix(counts_df[, all_samples]))

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = col_data,
    design    = ~ condition
  )
  dds <- DESeq(dds)

  res <- results(
    dds,
    contrast = c("condition", condition_a, condition_b)
  )

} else if (mode == "timeseries") {

  sample_columns <- unlist(p$sample_columns)
  factors        <- p$factors           # named list factor -> [per-sample vals]
  formula_full   <- as.formula(p$formula)
  formula_reduced <- as.formula(p$reduced)

  col_data <- as.data.frame(
    lapply(factors, function(vals) factor(unlist(vals))),
    row.names = sample_columns
  )
  count_mat <- round(as.matrix(counts_df[, sample_columns]))

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = col_data,
    design    = formula_full
  )
  dds <- DESeq(dds, test = "LRT", reduced = formula_reduced)
  res <- results(dds)

} else {
  stop(paste("Unknown mode:", mode), call. = FALSE)
}

# ---------------------------------------------------------------------------
# 4.  Rename columns and write output
# ---------------------------------------------------------------------------
df_res <- as.data.frame(res)
df_res <- df_res[rownames(counts_df), ]           # restore original row order

for (r_col in names(col_map)) {
  pipeline_col <- col_map[[r_col]]
  if (r_col %in% colnames(df_res)) {
    colnames(df_res)[colnames(df_res) == r_col] <- pipeline_col
  }
}

# Ensure output directories exist
dir.create(dirname(output_tsv),     recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_parquet), recursive = TRUE, showWarnings = FALSE)

message("Writing TSV     -> ", output_tsv)
write.table(df_res, file = output_tsv, sep = "\t", quote = FALSE, col.names = NA)

message("Writing Parquet -> ", output_parquet)
write_parquet(df_res, output_parquet)

message("DESeq2 done.")