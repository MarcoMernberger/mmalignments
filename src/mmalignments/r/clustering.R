#!/usr/bin/env Rscript
# clustering.R — Unsupervised clustering, PCA, scree/elbow, heatmap, volcano plots
#
# Usage:  Rscript clustering.R <params.json>
#
# The params JSON must contain:
#   mode : "pca" | "scree" | "heatmap" | "volcano"
#
# ---- mode == "pca" ----
#   counts_parquet   : path to count matrix (Parquet, genes x samples)
#   metadata_csv     : (optional) path to CSV with sample metadata
#   color_by         : (optional) metadata column for point colour
#   shape_by         : (optional) metadata column for point shape
#   show_labels      : bool — show labels on initial HTML render
#   label_spacing    : numeric — ggrepel box.padding for static PDF
#   pcs              : list of [pc_x, pc_y] pairs, e.g. [[1,2],[1,3]]
#   n_components     : int — number of PCs to compute
#   scale            : bool — scale genes to unit variance
#   center           : bool — center genes
#   comparison_name  : label for titles / file names
#   output_html      : path for interactive plotly HTML
#   output_pdf       : path for static ggplot2+ggrepel PDF
#   output_parquet   : path for PCA coords parquet
#
# ---- mode == "scree" ----
#   counts_parquet   : path to count matrix (Parquet)
#   n_components     : int — number of PCs to show
#   scale, center    : bool
#   comparison_name  : label
#   output_pdf       : path for PDF (scree + cumulative variance)
#
# ---- mode == "heatmap" ----
#   counts_parquet             : path to count matrix (Parquet)
#   metadata_csv               : (optional) path to metadata CSV
#   color_by                   : (optional) character vector of metadata columns
#   top_n_genes                : int — most-variable genes to show
#   clustering_distance_rows   : pheatmap distance metric for rows
#   clustering_distance_cols   : pheatmap distance metric for cols
#   clustering_method          : pheatmap linkage method
#   scale_rows                 : bool — z-score rows
#   comparison_name            : label
#   output_pdf                 : path for PDF
#
# ---- mode == "volcano" ----
#   results_parquet  : path to differential expression Parquet
#   logfc_col        : column name for log2FC
#   fdr_col          : column name for FDR / adjusted p-value
#   p_col            : (optional) column name for raw p-value (y-axis)
#   logfc_threshold  : numeric threshold for |log2FC|
#   fdr_threshold    : numeric threshold for FDR
#   label_top_n      : int — number of top genes to label
#   label_col        : (optional) column to use as gene labels
#   label_spacing    : numeric — ggrepel box.padding
#   show_labels      : bool — show labels on initial HTML render
#   comparison_name  : label
#   output_html      : path for interactive plotly HTML
#   output_pdf       : path for static PDF

suppressPackageStartupMessages({
  library(arrow)
  library(jsonlite)
  library(ggplot2)
  library(ggrepel)
  library(plotly)
  library(htmlwidgets)
  library(gridExtra)
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

read_counts <- function(counts_parquet) {
  df <- as.data.frame(read_parquet(counts_parquet))
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
  df
}

read_results <- function(results_parquet) {
  df <- as.data.frame(read_parquet(results_parquet))
  if (!is.null(df[[1]]) && !is.numeric(df[[1]])) {
    rownames(df) <- df[[1]]
    df[[1]] <- NULL
  }
  df
}

read_metadata <- function(metadata_csv, sample_names) {
  if (is.null(metadata_csv) || is.na(metadata_csv)) return(NULL)
  meta <- read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE)
  # first column must contain sample names
  rownames(meta) <- meta[[1]]
  meta[[1]] <- NULL
  meta[sample_names, , drop = FALSE]
}

ensure_dir <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
}

# Build a discrete colour palette that is colour-blind friendly
cb_palette <- function(n) {
  base <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#A65628", "#F781BF", "#999999",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"
  )
  if (n <= length(base)) return(base[seq_len(n)])
  colorRampPalette(base)(n)
}

# Estimate elbow via max distance from diagonal (kneedle approximation)
find_elbow <- function(y) {
  n <- length(y)
  if (n < 3) return(1L)
  straight <- seq(y[1], y[n], length.out = n)
  which.max(y - straight)
}

# ---------------------------------------------------------------------------
# Mode: PCA
# ---------------------------------------------------------------------------

run_pca <- function(p) {
  message("=== PCA ===")
  counts_df <- read_counts(p$counts_parquet)

  # Transpose: samples x genes for prcomp
  mat <- t(as.matrix(counts_df))

  # Remove zero-variance genes
  gene_var <- apply(mat, 2, var)
  mat <- mat[, gene_var > 0, drop = FALSE]

  n_comp <- min(p$n_components, nrow(mat) - 1, ncol(mat))
  message(sprintf("Running prcomp on %d samples x %d genes, n_comp=%d",
                  nrow(mat), ncol(mat), n_comp))

  pca_res <- prcomp(mat, scale. = isTRUE(p$scale), center = isTRUE(p$center),
                    rank. = n_comp)

  var_pct <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
  coords  <- as.data.frame(pca_res$x)
  coords$sample <- rownames(coords)

  # Merge metadata
  meta <- read_metadata(p$metadata_csv, rownames(mat))
  pca_df <- coords
  if (!is.null(meta)) {
    common <- intersect(rownames(meta), pca_df$sample)
    pca_df <- merge(pca_df, meta[common, , drop = FALSE],
                    by.x = "sample", by.y = "row.names", all.x = TRUE)
  }

  # Resolve color and shape columns
  color_col <- if (!is.null(p$color_by) && p$color_by %in% colnames(pca_df)) {
    p$color_by
  } else {
    NULL
  }
  shape_col <- if (!is.null(p$shape_by) && p$shape_by %in% colnames(pca_df)) {
    p$shape_by
  } else {
    NULL
  }

  pcs_list <- p$pcs  # list of [pc_x, pc_y]

  # ------------------------------------------------------------------
  # Static PDF (ggplot2 + ggrepel)
  # ------------------------------------------------------------------
  ensure_dir(p$output_pdf)
  ggplots <- lapply(pcs_list, function(pair) {
    pcx <- pair[[1]]; pcy <- pair[[2]]
    xcol <- paste0("PC", pcx); ycol <- paste0("PC", pcy)
    if (!xcol %in% colnames(pca_df) || !ycol %in% colnames(pca_df)) return(NULL)

    aes_base <- aes_string(x = xcol, y = ycol)
    g <- ggplot(pca_df, aes_base) +
      theme_bw(base_size = 11) +
      labs(
        title = p$comparison_name,
        x = sprintf("PC%d (%.1f%%)", pcx, var_pct[pcx]),
        y = sprintf("PC%d (%.1f%%)", pcy, var_pct[pcy])
      )

    if (!is.null(color_col) && !is.null(shape_col)) {
      g <- g + geom_point(aes_string(color = color_col, shape = shape_col),
                           size = 3, alpha = 0.85)
    } else if (!is.null(color_col)) {
      g <- g + geom_point(aes_string(color = color_col), size = 3, alpha = 0.85)
    } else if (!is.null(shape_col)) {
      g <- g + geom_point(aes_string(shape = shape_col), size = 3, alpha = 0.85)
    } else {
      g <- g + geom_point(size = 3, alpha = 0.85)
    }

    label_aes <- if (!is.null(color_col)) {
      aes_string(label = "sample", color = color_col)
    } else {
      aes_string(label = "sample")
    }
    g <- g + geom_text_repel(
      label_aes,
      size          = 2.8,
      box.padding   = p$label_spacing,
      max.overlaps  = Inf,
      show.legend   = FALSE
    )
    g
  })
  ggplots <- Filter(Negate(is.null), ggplots)

  ncols_pdf <- min(length(ggplots), 2)
  nrows_pdf <- ceiling(length(ggplots) / ncols_pdf)
  pdf(p$output_pdf, width = 7 * ncols_pdf, height = 6 * nrows_pdf)
  do.call(grid.arrange, c(ggplots, list(ncol = ncols_pdf)))
  dev.off()
  message("Wrote PDF: ", p$output_pdf)

  # ------------------------------------------------------------------
  # Interactive HTML (plotly)
  # ------------------------------------------------------------------
  ensure_dir(p$output_html)

  subplot_list <- lapply(pcs_list, function(pair) {
    pcx <- pair[[1]]; pcy <- pair[[2]]
    xcol <- paste0("PC", pcx); ycol <- paste0("PC", pcy)
    if (!xcol %in% colnames(pca_df) || !ycol %in% colnames(pca_df)) return(NULL)

    # Hover text
    hover_lines <- paste0("Sample: ", pca_df$sample)
    if (!is.null(color_col)) {
      hover_lines <- paste0(hover_lines, "<br>", color_col, ": ", pca_df[[color_col]])
    }
    hover_lines <- paste0(
      hover_lines,
      "<br>", sprintf("PC%d: %.3f", pcx, pca_df[[xcol]]),
      "<br>", sprintf("PC%d: %.3f", pcy, pca_df[[ycol]])
    )

    color_vals <- if (!is.null(color_col)) pca_df[[color_col]] else "Samples"

    # --- Trace 0: markers ---
    trace_markers <- list(
      x          = pca_df[[xcol]],
      y          = pca_df[[ycol]],
      text       = hover_lines,
      hoverinfo  = "text",
      mode       = "markers",
      type       = "scatter",
      name       = "Samples",
      showlegend = TRUE,
      marker     = list(size = 9, opacity = 0.85)
    )
    if (!is.null(color_col)) {
      lvls <- unique(pca_df[[color_col]])
      pal  <- cb_palette(length(lvls))
      col_map <- setNames(pal, lvls)
      trace_markers$marker$color <- col_map[as.character(pca_df[[color_col]])]
    }

    # --- Trace 1: text labels (separate trace for independent font-size slider) ---
    default_size <- if (isTRUE(p$show_labels)) 11 else 0
    trace_labels <- list(
      x          = pca_df[[xcol]],
      y          = pca_df[[ycol]],
      text       = pca_df$sample,
      mode       = "text",
      type       = "scatter",
      hoverinfo  = "skip",
      showlegend = FALSE,
      name       = "Labels",
      textfont   = list(size = default_size, color = "rgba(50,50,50,0.9)"),
      textposition = "top center"
    )

    xlab <- sprintf("PC%d (%.1f%%)", pcx, var_pct[pcx])
    ylab <- sprintf("PC%d (%.1f%%)", pcy, var_pct[pcy])

    list(
      traces = list(trace_markers, trace_labels),
      xlab   = xlab,
      ylab   = ylab,
      default_size = default_size
    )
  })
  subplot_list <- Filter(Negate(is.null), subplot_list)

  # Font-size slider steps (applied to ALL label traces via restyle)
  label_sizes <- c(0, 7, 9, 11, 13, 16)
  default_size_all <- subplot_list[[1]]$default_size

  # Index of label traces among all traces:
  # subplot i contributes traces at positions (2i-2, 2i-1) => label trace at 2i-1
  label_trace_indices <- seq(1, 2 * length(subplot_list) - 1, by = 2)  # 0-based

  slider_steps <- lapply(label_sizes, function(sz) {
    list(
      method = "restyle",
      args   = list("textfont.size", sz, as.list(label_trace_indices)),
      label  = if (sz == 0) "Off" else as.character(sz)
    )
  })
  active_step <- which(label_sizes == default_size_all) - 1
  if (length(active_step) == 0) active_step <- 3L

  sliders_layout <- list(list(
    active       = active_step,
    currentvalue = list(prefix = "Label size: ", font = list(size = 13)),
    pad          = list(t = 50),
    x            = 0.0,
    len          = 0.45,
    steps        = slider_steps
  ))

  # Build figure
  fig <- plot_ly()
  for (sp in subplot_list) {
    for (tr in sp$traces) {
      fig <- do.call(add_trace, c(list(fig), tr))
    }
  }

  # Use subplot if more than one PC pair
  if (length(subplot_list) > 1) {
    # Rebuild with proper subplot layout
    sub_figs <- lapply(subplot_list, function(sp) {
      f <- plot_ly()
      for (tr in sp$traces) {
        f <- do.call(add_trace, c(list(f), tr))
      }
      f <- layout(f, xaxis = list(title = sp$xlab),
                     yaxis = list(title = sp$ylab))
      f
    })
    ncols_html <- min(length(sub_figs), 2)
    fig <- subplot(sub_figs, nrows = ceiling(length(sub_figs) / ncols_html),
                   shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE)
  } else {
    fig <- layout(fig,
      xaxis = list(title = subplot_list[[1]]$xlab),
      yaxis = list(title = subplot_list[[1]]$ylab)
    )
  }

  fig <- layout(fig,
    title   = list(text = p$comparison_name, x = 0.5),
    sliders = sliders_layout,
    margin  = list(b = 120)
  )

  saveWidget(fig, file = p$output_html, selfcontained = TRUE,
             title = paste("PCA —", p$comparison_name))
  message("Wrote HTML: ", p$output_html)

  # ------------------------------------------------------------------
  # Parquet: PCA coordinates
  # ------------------------------------------------------------------
  ensure_dir(p$output_parquet)
  write_parquet(coords, p$output_parquet)
  message("Wrote Parquet: ", p$output_parquet)
}

# ---------------------------------------------------------------------------
# Mode: scree
# ---------------------------------------------------------------------------

run_scree <- function(p) {
  message("=== Scree / Elbow ===")
  counts_df <- read_counts(p$counts_parquet)

  mat      <- t(as.matrix(counts_df))
  gene_var <- apply(mat, 2, var)
  mat      <- mat[, gene_var > 0, drop = FALSE]

  n_comp <- min(p$n_components, nrow(mat) - 1, ncol(mat))
  pca_res <- prcomp(mat, scale. = isTRUE(p$scale), center = isTRUE(p$center),
                    rank. = n_comp)

  n_show  <- min(p$n_components, length(pca_res$sdev))
  var_pct <- (pca_res$sdev[seq_len(n_show)]^2 / sum(pca_res$sdev^2)) * 100
  cum_pct <- cumsum(pca_res$sdev[seq_len(n_show)]^2) / sum(pca_res$sdev^2) * 100

  var_df <- data.frame(PC = seq_len(n_show), variance = var_pct, cumulative = cum_pct)
  elbow  <- find_elbow(var_pct)

  title_str <- p$comparison_name

  p_scree <- ggplot(var_df, aes(x = PC, y = variance)) +
    geom_col(fill = "#377EB8", colour = "white", width = 0.8) +
    geom_vline(xintercept = elbow, colour = "#E41A1C", linetype = "dashed",
               linewidth  = 0.8) +
    annotate("text", x = elbow + 0.5, y = max(var_pct) * 0.92,
             label = sprintf("Elbow\n(PC%d)", elbow),
             hjust = 0, colour = "#E41A1C", size = 3.5) +
    scale_x_continuous(breaks = seq_len(n_show)) +
    theme_bw(base_size = 11) +
    labs(title = paste("Scree —", title_str),
         x = "Principal Component",
         y = "Variance explained (%)")

  p_cum <- ggplot(var_df, aes(x = PC, y = cumulative)) +
    geom_line(colour = "#377EB8", linewidth = 0.9) +
    geom_point(colour = "#377EB8", size = 2.5) +
    geom_vline(xintercept = elbow, colour = "#E41A1C", linetype = "dashed",
               linewidth  = 0.8) +
    geom_hline(yintercept = 80, colour = "gray50", linetype = "dotted",
               linewidth  = 0.7) +
    annotate("text", x = max(n_show) * 0.65, y = 81.5,
             label = "80 %", colour = "gray40", size = 3.2) +
    scale_x_continuous(breaks = seq_len(n_show)) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_bw(base_size = 11) +
    labs(title = paste("Cumulative variance —", title_str),
         x = "Principal Component",
         y = "Cumulative variance (%)")

  ensure_dir(p$output_pdf)
  pdf(p$output_pdf, width = 14, height = 5.5)
  grid.arrange(p_scree, p_cum, ncol = 2)
  dev.off()
  message("Wrote PDF: ", p$output_pdf)
}

# ---------------------------------------------------------------------------
# Mode: heatmap
# ---------------------------------------------------------------------------

run_heatmap <- function(p) {
  message("=== Heatmap ===")
  suppressPackageStartupMessages(library(pheatmap))

  counts_df <- read_counts(p$counts_parquet)
  mat       <- as.matrix(counts_df)

  # Select top variable genes by IQR
  top_n   <- min(p$top_n_genes, nrow(mat))
  iqr_val <- apply(mat, 1, IQR)
  top_genes <- names(sort(iqr_val, decreasing = TRUE))[seq_len(top_n)]
  sub_mat <- mat[top_genes, , drop = FALSE]

  # Z-score scale rows
  if (isTRUE(p$scale_rows)) {
    sub_mat <- t(scale(t(sub_mat)))
    sub_mat[is.nan(sub_mat)] <- 0
    sub_mat[is.infinite(sub_mat)] <- 0
  }

  # Sample annotation
  annotation_col <- NA
  meta <- read_metadata(p$metadata_csv, colnames(sub_mat))
  if (!is.null(meta) && !is.null(p$color_by) && length(p$color_by) > 0) {
    cols_present <- intersect(p$color_by, colnames(meta))
    if (length(cols_present) > 0) {
      annotation_col <- meta[colnames(sub_mat), cols_present, drop = FALSE]
    }
  }

  show_rownames <- (top_n <= 50)
  fontsize_row  <- if (top_n <= 50) 8 else 6
  height_in     <- max(8, min(24, top_n / 20 + 4))

  ensure_dir(p$output_pdf)
  pheatmap(
    sub_mat,
    annotation_col          = annotation_col,
    clustering_distance_rows = p$clustering_distance_rows,
    clustering_distance_cols = p$clustering_distance_cols,
    clustering_method        = p$clustering_method,
    show_rownames            = show_rownames,
    fontsize                 = 9,
    fontsize_row             = fontsize_row,
    main                     = p$comparison_name,
    filename                 = p$output_pdf,
    width                    = 12,
    height                   = height_in
  )
  message("Wrote PDF: ", p$output_pdf)
}

# ---------------------------------------------------------------------------
# Mode: volcano
# ---------------------------------------------------------------------------

run_volcano <- function(p) {
  message("=== Volcano ===")

  res_df <- read_results(p$results_parquet)

  logfc_col <- p$logfc_col
  fdr_col   <- p$fdr_col
  y_col     <- if (!is.null(p$p_col) && p$p_col %in% colnames(res_df)) p$p_col else fdr_col

  for (col in c(logfc_col, fdr_col)) {
    if (!col %in% colnames(res_df)) {
      stop(sprintf("Column '%s' not found in results table. Available: %s",
                   col, paste(colnames(res_df), collapse = ", ")))
    }
  }

  lfc_thr <- p$logfc_threshold
  fdr_thr <- p$fdr_threshold

  res_df$sig <- "n.s."
  res_df$sig[!is.na(res_df[[fdr_col]]) &
               res_df[[fdr_col]]  < fdr_thr &
               res_df[[logfc_col]] >  lfc_thr] <- "up"
  res_df$sig[!is.na(res_df[[fdr_col]]) &
               res_df[[fdr_col]]  < fdr_thr &
               res_df[[logfc_col]] < -lfc_thr] <- "down"

  # −log10 y-axis
  res_df$.y_val <- -log10(pmax(res_df[[y_col]], .Machine$double.xmin))

  # Gene labels
  label_col <- if (!is.null(p$label_col) && p$label_col %in% colnames(res_df)) {
    p$label_col
  } else {
    ".gene_label"
  }
  res_df$.gene_label <- rownames(res_df)

  # Top genes to label: rank by FDR * 1/(|logFC|+1e-9)
  score_rank <- !is.na(res_df[[fdr_col]]) & res_df[[fdr_col]] > 0
  res_df$.rank_score <- Inf
  res_df$.rank_score[score_rank] <- (
    res_df[[fdr_col]][score_rank] / (abs(res_df[[logfc_col]][score_rank]) + 1e-9)
  )
  top_idx    <- order(res_df$.rank_score)[seq_len(min(p$label_top_n, nrow(res_df)))]
  label_df   <- res_df[top_idx, ]

  sig_colours <- c("up" = "#D73027", "down" = "#4575B4", "n.s." = "#AAAAAA")
  y_axis_label <- sprintf("-log10(%s)", y_col)

  # ------------------------------------------------------------------
  # Static PDF (ggplot2 + ggrepel)
  # ------------------------------------------------------------------
  ensure_dir(p$output_pdf)

  g <- ggplot(res_df, aes_string(x = logfc_col, y = ".y_val", colour = "sig")) +
    geom_point(alpha = 0.55, size = 1.4) +
    scale_colour_manual(values = sig_colours, name = "Significance") +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr),
               linetype = "dashed", colour = "gray40", linewidth = 0.6) +
    geom_hline(yintercept = -log10(fdr_thr),
               linetype = "dashed", colour = "gray40", linewidth = 0.6) +
    geom_text_repel(
      data          = label_df,
      aes_string(x  = logfc_col, y = ".y_val", label = label_col),
      colour        = "black",
      size          = 2.7,
      box.padding   = p$label_spacing,
      max.overlaps  = Inf,
      segment.alpha = 0.5,
      inherit.aes   = FALSE
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = p$comparison_name,
      x     = sprintf("log\u2082 fold change  (%s)", logfc_col),
      y     = y_axis_label
    )

  ggsave(p$output_pdf, g, width = 8, height = 7, device = "pdf")
  message("Wrote PDF: ", p$output_pdf)

  # ------------------------------------------------------------------
  # Interactive HTML (plotly)
  # ------------------------------------------------------------------
  ensure_dir(p$output_html)

  res_df$hover_text <- paste0(
    "Gene: ", res_df$.gene_label,
    "<br>log2FC: ", round(res_df[[logfc_col]], 3),
    "<br>", y_col, ": ", signif(res_df[[y_col]], 3),
    "<br>Significance: ", res_df$sig
  )

  default_size <- if (isTRUE(p$show_labels)) 11 else 0

  fig <- plot_ly()

  for (grp in c("up", "down", "n.s.")) {
    sub <- res_df[res_df$sig == grp, ]
    if (nrow(sub) == 0) next
    fig <- add_trace(fig,
      x         = sub[[logfc_col]],
      y         = sub$.y_val,
      text      = sub$hover_text,
      hoverinfo = "text",
      mode      = "markers",
      type      = "scatter",
      name      = grp,
      marker    = list(color = sig_colours[[grp]], size = 7, opacity = 0.6)
    )
  }

  # Label trace (index = number of sig groups so far, 0-based)
  n_group_traces <- length(unique(res_df$sig[res_df$sig %in% c("up", "down", "n.s.")]))
  label_trace_idx <- n_group_traces  # 0-based index of the labels trace

  fig <- add_trace(fig,
    x          = label_df[[logfc_col]],
    y          = label_df$.y_val,
    text       = label_df[[label_col]],
    mode       = "text",
    type       = "scatter",
    hoverinfo  = "skip",
    showlegend = FALSE,
    name       = "Labels",
    textfont   = list(size = default_size, color = "rgba(30,30,30,0.85)"),
    textposition = "top center"
  )

  # Font-size slider for labels
  label_sizes  <- c(0, 7, 9, 11, 13, 16)
  active_step  <- which(label_sizes == default_size) - 1
  if (length(active_step) == 0) active_step <- 3L

  slider_steps <- lapply(label_sizes, function(sz) {
    list(
      method = "restyle",
      args   = list("textfont.size", sz, list(label_trace_idx)),
      label  = if (sz == 0) "Off" else as.character(sz)
    )
  })

  fig <- layout(fig,
    title   = list(text = p$comparison_name, x = 0.5),
    xaxis   = list(title = sprintf("log\u2082 fold change  (%s)", logfc_col),
                   zeroline = TRUE, zerolinecolor = "#cccccc"),
    yaxis   = list(title = y_axis_label),
    shapes  = list(
      list(type = "line", x0 = -lfc_thr, x1 = -lfc_thr,
           y0 = 0, y1 = 1, yref = "paper",
           line = list(dash = "dash", color = "gray", width = 1)),
      list(type = "line", x0 = lfc_thr, x1 = lfc_thr,
           y0 = 0, y1 = 1, yref = "paper",
           line = list(dash = "dash", color = "gray", width = 1)),
      list(type = "line", x0 = 0, x1 = 1, xref = "paper",
           y0 = -log10(fdr_thr), y1 = -log10(fdr_thr),
           line = list(dash = "dash", color = "gray", width = 1))
    ),
    sliders = list(list(
      active       = active_step,
      currentvalue = list(prefix = "Label size: ", font = list(size = 13)),
      pad          = list(t = 50),
      x            = 0.0,
      len          = 0.45,
      steps        = slider_steps
    )),
    margin  = list(b = 120)
  )

  saveWidget(fig, file = p$output_html, selfcontained = TRUE,
             title = paste("Volcano —", p$comparison_name))
  message("Wrote HTML: ", p$output_html)
}

# ---------------------------------------------------------------------------
# Main dispatcher
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript clustering.R <params.json>", call. = FALSE)
}

p <- fromJSON(args[1], simplifyVector = FALSE)

switch(
  p$mode,
  pca     = run_pca(p),
  scree   = run_scree(p),
  heatmap = run_heatmap(p),
  volcano = run_volcano(p),
  stop(sprintf("Unknown mode '%s'. Use: pca, scree, heatmap, volcano", p$mode),
       call. = FALSE)
)

message("Done.")
