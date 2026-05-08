#!/usr/bin/env Rscript
# edger_unpaired.R — edgeR exactTest for two-group unpaired comparisons
#
# Usage:  Rscript edger_unpaired.R <params.json>
#
# Expected params JSON keys:
#   counts_parquet          : path to count matrix (Parquet)
#   condition_a             : treatment condition label
#   condition_b             : reference condition label
#   columns_a               : list of sample column names for condition_a
#   columns_b               : list of sample column names for condition_b
#   comparison_name         : human-readable label (for logging)
#   library_sizes           : list of floats, or null (use column sums)
#   manual_dispersion_value : float, fallback dispersion for n=1 replicates
#   column_map              : named object R-col -> pipeline-col
#   output_tsv              : output path
#   output_parquet          : output path

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
  stop("Usage: Rscript edger_unpaired.R <params.json>", call. = FALSE)
}

p <- fromJSON(args[1], simplifyVector = FALSE)

columns_a   <- unlist(p$columns_a)
columns_b   <- unlist(p$columns_b)
col_map     <- p$column_map
output_tsv     <- p$output_tsv
output_parquet <- p$output_parquet

# ---------------------------------------------------------------------------
# 2.  Load counts
# ---------------------------------------------------------------------------
message("Reading counts from ", p$counts_parquet)
counts_df <- as.data.frame(read_parquet(p$counts_parquet))
rownames(counts_df) <- counts_df[[1]]
counts_df[[1]] <- NULL

all_cols  <- c(columns_a, columns_b)
count_mat <- round(as.matrix(counts_df[, all_cols]))

# Rename to X_0, X_1, … to avoid edgeR complaints about special chars
renamed <- paste0("X_", seq_len(ncol(count_mat)) - 1)
colnames(count_mat) <- renamed

# ---------------------------------------------------------------------------
# 3.  Build DGEList
# ---------------------------------------------------------------------------
n_a <- length(columns_a)
n_b <- length(columns_b)

lib_sizes <- if (!is.null(p$library_sizes)) {
  unlist(p$library_sizes)
} else {
  colSums(count_mat)
}

samples_df <- data.frame(
  group    = factor(c(rep("z", n_a), rep("x", n_b))),  # z=treatment, x=ref
  lib.size = lib_sizes,
  row.names = renamed
)

dge <- DGEList(counts = count_mat, samples = samples_df)
dge <- calcNormFactors(dge)

# ---------------------------------------------------------------------------
# 4.  Dispersion + exactTest
# ---------------------------------------------------------------------------
if (n_a == 1 && n_b == 1) {
  message(
    "Single replicate per condition — using manual dispersion ",
    p$manual_dispersion_value,
    ".\n",
    "Consider adding more replicates or using logFC only."
  )
  exact_tested <- exactTest(
    dge,
    dispersion = (p$manual_dispersion_value)^2
  )
} else {
  disp        <- estimateDisp(dge, robust = TRUE)
  exact_tested <- exactTest(disp)
}

res <- topTags(exact_tested, n = nrow(count_mat), sort.by = "none")
df_res <- as.data.frame(res$table)
df_res <- df_res[rownames(counts_df), ]     # restore original gene order

# ---------------------------------------------------------------------------
# 5.  Rename columns
# ---------------------------------------------------------------------------
for (r_col in names(col_map)) {
  pipeline_col <- col_map[[r_col]]
  if (r_col %in% colnames(df_res)) {
    colnames(df_res)[colnames(df_res) == r_col] <- pipeline_col
  }
}

# ---------------------------------------------------------------------------
# 6.  Write output
# ---------------------------------------------------------------------------
dir.create(dirname(output_tsv),     recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_parquet), recursive = TRUE, showWarnings = FALSE)

message("Writing TSV     -> ", output_tsv)
write.table(df_res, file = output_tsv, sep = "\t", quote = FALSE, col.names = NA)

message("Writing Parquet -> ", output_parquet)
write_parquet(df_res, output_parquet)

message("edgeR done.")