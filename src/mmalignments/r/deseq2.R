#!/usr/bin/env Rscript
# deseq2.R — DESeq2 differential expression and normalisation functions
#
# Called in-process via rpy2 (RScriptInternal).  Each exported function
# accepts R data.frames / vectors and returns an R data.frame.
# File I/O (TSV / Parquet) is handled entirely on the Python side.

library(DESeq2)

# ---------------------------------------------------------------------------
# Helper: detect sample columns (numeric columns of a data.frame)
# ---------------------------------------------------------------------------
.detect_sample_columns <- function(df) {
    colnames(df)[sapply(df, is.numeric)]
}

# ---------------------------------------------------------------------------
# counts_to_vst
#   Variance-stabilising transformation (blind or design-aware).
#
# Parameters
#   counts_df      : data.frame  genes × (annotations + samples)
#   sample_columns : character vector of count column names, or NULL (→ auto)
#   blind          : logical, use blind dispersion estimation (default TRUE)
#   fitType        : "parametric" | "local" | "mean"
#
# Returns
#   data.frame  genes × sample_columns  with VST values
# ---------------------------------------------------------------------------
counts_to_vst <- function(
    counts_df,
    sample_columns = NULL,
    blind          = TRUE,
    fitType        = "parametric"
) {
    df <- counts_df

    if (is.null(sample_columns)) {
        sample_columns <- .detect_sample_columns(df)
    }

    count_mat         <- round(as.matrix(df[, sample_columns, drop = FALSE]))
    rownames(count_mat) <- rownames(df)

    col_data <- data.frame(
        row.names = sample_columns,
        condition = factor(rep("all", length(sample_columns)))
    )

    dds <- DESeqDataSetFromMatrix(
        countData = count_mat,
        colData   = col_data,
        design    = ~1
    )

    vsd    <- vst(dds, blind = blind, fitType = fitType)
    vst_df <- as.data.frame(assay(vsd))

    return(vst_df)
}

# ---------------------------------------------------------------------------
# counts_to_rlog
#   Regularised log transformation.
#
# Parameters
#   counts_df      : data.frame  genes × (annotations + samples)
#   sample_columns : character vector of count column names, or NULL (→ auto)
#   blind          : logical (default TRUE)
#
# Returns
#   data.frame  genes × sample_columns  with rlog values
# ---------------------------------------------------------------------------
counts_to_rlog <- function(
    counts_df,
    sample_columns = NULL,
    blind          = TRUE
) {
    df <- counts_df

    if (is.null(sample_columns)) {
        sample_columns <- .detect_sample_columns(df)
    }

    count_mat         <- round(as.matrix(df[, sample_columns, drop = FALSE]))
    rownames(count_mat) <- rownames(df)

    col_data <- data.frame(
        row.names = sample_columns,
        condition = factor(rep("all", length(sample_columns)))
    )

    dds     <- DESeqDataSetFromMatrix(
        countData = count_mat,
        colData   = col_data,
        design    = ~1
    )

    rld      <- rlog(dds, blind = blind)
    rlog_df  <- as.data.frame(assay(rld))

    return(rlog_df)
}

# ---------------------------------------------------------------------------
# deseq2_unpaired
#   Two-group Wald-test comparison.
#
# Parameters
#   counts_df        : data.frame  genes × (annotations + samples)
#   condition_a      : treatment condition label (numerator)
#   condition_b      : reference condition label (denominator)
#   model_conditions : named list  condition -> character vector of sample cols
#   column_map       : named list  DESeq2 col -> pipeline col (for renaming)
#
# Returns
#   data.frame with per-gene statistics (column names from DESeq2, or renamed
#   if column_map is provided)
# ---------------------------------------------------------------------------
deseq2_unpaired <- function(
    counts_df,
    condition_a,
    condition_b,
    model_conditions,
    column_map = list()
) {
    cat("=== counts_df ===\n")
    print(head(counts_df))
    cat("\nDimensions:", nrow(counts_df), "genes x", ncol(counts_df), "samples\n")
    cat("Columns:\n")
    print(colnames(counts_df))
    cat("First gene IDs:\n")
    print(head(rownames(counts_df)))

    # build sample metadata
    all_samples <- c()
    all_labels  <- c()

    for (cond in names(model_conditions)) {
        cols <- as.character(unlist(model_conditions[[cond]]))
        all_samples <- c(all_samples, cols)
        all_labels  <- c(all_labels, rep(cond, length(cols)))
    }

    cat("\n=== model_conditions ===\n")
    print(model_conditions)

    cat("\n=== all_samples ===\n")
    print(all_samples)

    cat("\n=== missing samples ===\n")
    print(setdiff(all_samples, colnames(counts_df)))

    extra_levels <- setdiff(
        names(model_conditions),
        c(condition_a, condition_b)
    )

    col_data <- data.frame(
        condition = factor(
            all_labels,
            levels = c(condition_a, condition_b, extra_levels)
        ),
        row.names = all_samples
    )

    count_mat <- round(as.matrix(counts_df[, all_samples, drop = FALSE]))

    dds <- DESeqDataSetFromMatrix(
        countData = count_mat,
        colData   = col_data,
        design    = ~condition
    )

    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", condition_a, condition_b))
    df_res <- as.data.frame(res)
    df_res <- df_res[rownames(counts_df), ]

    for (r_col in names(column_map)) {
        pipeline_col <- as.character(column_map[[r_col]])
        if (r_col %in% colnames(df_res)) {
            colnames(df_res)[colnames(df_res) == r_col] <- pipeline_col
        }
    }
    sf <- sizeFactors(dds)
    result_bundle <- list(
        results = as.data.frame(res),
        size_factors = as.data.frame(
            sample = names(sizeFactors(dds)),
            size_factor = as.numeric(sizeFactors(dds))
        )
    )
    # df_res
    return(result_bundle)
}

# ---------------------------------------------------------------------------
# deseq2_timeseries
#   Multi-factor / time-series LRT comparison.
#
# Parameters
#   counts_df      : data.frame  genes × (annotations + samples)
#   sample_columns : character vector of sample column names (ordered)
#   factors        : named list  factor_name -> per-sample values
#   formula        : full model formula string, e.g. "~ condition + time"
#   reduced        : reduced model formula string, e.g. "~ condition"
#   column_map     : named list  DESeq2 col -> pipeline col
#
# Returns
#   data.frame with per-gene LRT statistics
# ---------------------------------------------------------------------------
deseq2_timeseries <- function(
    counts_df,
    sample_columns,
    factors,
    formula,
    reduced,
    column_map = list()
) {
    sample_columns <- as.character(unlist(sample_columns))

    col_data <- as.data.frame(
        lapply(factors, function(vals) factor(as.character(unlist(vals)))),
        row.names = sample_columns
    )

    count_mat         <- round(as.matrix(counts_df[, sample_columns, drop = FALSE]))
    rownames(count_mat) <- rownames(counts_df)

    dds <- DESeqDataSetFromMatrix(
        countData = count_mat,
        colData   = col_data,
        design    = as.formula(formula)
    )
    dds    <- DESeq(dds, test = "LRT", reduced = as.formula(reduced))
    res    <- results(dds)
    df_res <- as.data.frame(res)
    df_res <- df_res[rownames(counts_df), ]

    for (r_col in names(column_map)) {
        pipeline_col <- as.character(column_map[[r_col]])
        if (r_col %in% colnames(df_res)) {
            colnames(df_res)[colnames(df_res) == r_col] <- pipeline_col
        }
    }

    return(df_res)
}

