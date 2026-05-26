library("MAGeCKFlute")
library("ggplot2")


NormalizeMageckBetaTable <- function(
    df,
    gene_col    = "Gene",       # gene id column
    beta_pattern = "\\|beta$",  # how to find the beta columns
    raw_suffix = "_raw"         # suffix to rename unnormalized betas 
) {
    # package check
    if (!requireNamespace("MAGeCKFlute", quietly = TRUE)) {
        stop("Package 'MAGeCKFlute' is required but not installed.")
    }

    if (!gene_col %in% colnames(df)) {
        stop("gene_col '", gene_col, "' not found in input table.")
    }

    # get beta cols
    beta_cols <- grep(beta_pattern, colnames(df), value = TRUE)
    if (length(beta_cols) == 0) {
        stop("No beta columns found with pattern '", beta_pattern, "'.")
    }

    # beta frame for NormalizeBeta
    beta_df <- df[, c(gene_col, beta_cols), drop = FALSE]

    # normalize
    beta_norm <- MAGeCKFlute::NormalizeBeta(beta_df)
    beta_norm_only <- beta_norm[, beta_cols, drop = FALSE]

    # new colnames for original beta columns
    raw_colnames <- paste0(beta_cols, raw_suffix)
    colnames(df)[match(beta_cols, colnames(df))] <- raw_colnames    

    # merge
    df_out <- cbind(df, beta_norm_only)
    return(df_out)
}

RunMageckScatterView <- function(
    input_file,  # table with fold changes/beta scores
    data_type = c("mle", "rra"),  # "mle" = MLE table, "rra" = RRA / generic table
    output_dir,        # where to put the files, gets created
    filebase_name,     # filebase name for the output files
    x_col,
    y_col,
    gene_col = "id",   # MLE often "Gene", RRA often "id"
    sep = "\t",
    normalize = TRUE,  # for MLE: use NormalizeBeta, for RRA: not implemented yet
    top = 10,          # top x genes per group in ScatterView
    groups = c("bottomleft"),  # this controls selection! quadrants of interest for ScatterView & selection, NULL selects all
    # selection parameter based on the diagonal
    select = NULL,     # this controls selection! one of NULL, "positive", "negative", "both", "none" 
    neg_effect_cutoff = -0.4,
    pos_effect_cutoff = 0.4,
    delta_cutoff_k = 2,
    filter_fdr_x = FALSE,      # this controls selection! if True, select for significant genes in x_col
    filter_fdr_y = FALSE,      # this controls selection! if True, select for significant genes in y_col
    filter_groups = TRUE,
    fdr_cutoff = 0.05,     # FDR cutoff
    # plot parameters
    toplabels=TRUE,   # does not affect selection, only the plot
    label_selected_only=FALSE,
    xlab = NULL,
    ylab = NULL,
    jpeg_width = 20,
    jpeg_height = 15,
    # additional plot parameter for MLE
    auto_cut_diag = 2,
    auto_cut_x = 2,
    auto_cut_y = 2,
    # additional plot parameter for RRA
    x_cut = NULL, # c(-neg_effect_cutoff, pos_effect_cutoff),
    y_cut = NULL # c(-neg_effect_cutoff, pos_effect_cutoff)
) {

    infer_fdr_col <- function(colname) {
        if (grepl("\\|beta$", colname)) {
            # MLE
            return(sub("\\|beta$", "|wald-fdr", colname))
        }
        else if (grepl("\\|lfc$", colname)) {
            # RRA
            return(sub("\\|lfc$", "|fdr", colname))
        }
        else {
            # stop(sprintf("Could not infer FDR column from '%s'", colname))
            return(NULL)
        }
    }

    data_type <- match.arg(data_type)
    fdr_x_col <- infer_fdr_col(x_col)
    fdr_y_col <- infer_fdr_col(y_col)

    # required packages
    if (!requireNamespace("MAGeCKFlute", quietly = TRUE)) {
        stop("Package 'MAGeCKFlute' is required.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    # create output folder
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # read the data
    df <- utils::read.table(input_file, header = TRUE, sep = sep, check.names = FALSE)

    # check columns
    required_cols <- c(x_col, y_col, gene_col, fdr_x_col, fdr_y_col)
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
        stop(
            "Missing required columns in input data: ",
            paste(missing_cols, collapse = ", "),
            "  Present were: ",
            paste(colnames(df), collapse = ", ")
        )
    }

    # optional normalization, recommended for MLE
    df_norm <- df
    if (data_type == "mle" && normalize) {
        df_norm <- NormalizeMageckBetaTable(df)
    } else if (data_type == "rra" && normalize) {
        warning("normalize = TRUE is set for data_type='rra', but no normalization is implemented. Proceeding without changes.")
    }

    # ensure x/y are numeric
    for (col in c(x_col, y_col)) {
        if (!is.numeric(df_norm[[col]])) {
            suppressWarnings(
                df_norm[[col]] <- as.numeric(df_norm[[col]])
            )
        }
    }

    # default axis labels
    if (is.null(xlab)) xlab <- x_col
    if (is.null(ylab)) ylab <- y_col

    n <- nrow(df_norm)

    # lineare regression for RRA-delta
    formula_str <- as.formula(
        sprintf("`%s` ~ `%s`", y_col, x_col)
    )
    fit <- stats::lm(formula_str, data = df_norm)
    coefs <- stats::coef(fit)

    if (length(coefs) < 2 || all(is.na(coefs[-1]))) {
        stop("Could not estimate slope from lm(). Check that x_col and y_col are numeric and not all NA.")
    }

    slope_fit <- unname(coefs[2])

    # delta & Cutoffs
    if (data_type == "mle") {
        delta <- df_norm[[y_col]] - df_norm[[x_col]]
    } else {
        delta <- df_norm[[y_col]] - slope_fit * df_norm[[x_col]]
    }

    delta_cutoff <- MAGeCKFlute::CutoffCalling(delta, delta_cutoff_k)
    intercept_min <- -delta_cutoff
    intercept_max <-  delta_cutoff

    # diagonal classification (pos/neg/none)
    diagonal <- rep("none", n)
    diagonal[delta >  delta_cutoff] <- "pos"  # y stronger as x
    diagonal[delta < -delta_cutoff] <- "neg"  # y weaker as x
    diagonal <- factor(diagonal, levels = c("none", "pos", "neg"))

    # verbose description for diag
    diag_desc <- rep("no difference", n)
    diag_desc[delta >  delta_cutoff] <- paste0(y_col, " > ", x_col)
    diag_desc[delta < -delta_cutoff] <- paste0(y_col, " < ", x_col)

    # selection X and Y (pos/neg/none) - are the genes selected for in screen?
    sel_x_cat <- ifelse(df_norm[[x_col]] <= neg_effect_cutoff, "neg",
                   ifelse(df_norm[[x_col]] >= pos_effect_cutoff, "pos", "none"))
    sel_y_cat <- ifelse(df_norm[[y_col]] <= neg_effect_cutoff, "neg",
                   ifelse(df_norm[[y_col]] >= pos_effect_cutoff, "pos", "none"))

    # quadrants 9
    quadrant <- ifelse(sel_y_cat == "pos" & sel_x_cat == "neg",  "topleft",
                  ifelse(sel_y_cat == "pos" & sel_x_cat == "none", "topcenter",
                  ifelse(sel_y_cat == "pos" & sel_x_cat == "pos",  "topright",
                  ifelse(sel_y_cat == "none" & sel_x_cat == "neg",  "midleft",
                  ifelse(sel_y_cat == "none" & sel_x_cat == "none", "midcenter",
                  ifelse(sel_y_cat == "none" & sel_x_cat == "pos",  "midright",
                  ifelse(sel_y_cat == "neg" & sel_x_cat == "neg",  "bottomleft",
                  ifelse(sel_y_cat == "neg" & sel_x_cat == "none", "bottomcenter",
                  ifelse(sel_y_cat == "neg" & sel_x_cat == "pos",  "bottomright", "unknown")))))))))

    quadrant <- factor(quadrant,
                       levels = c("topleft","topcenter","topright",
                                  "midleft","midcenter","midright",
                                  "bottomleft","bottomcenter","bottomright",
                                  "unknown"))

    # FDR filtering
    sig_x <- rep(TRUE, n)
    sig_y <- rep(TRUE, n)

    if (!is.null(fdr_x_col)) {
        # FDR filter - if not significant_filter -> FALSE
        sig_x <- df_norm[[fdr_x_col]] <= fdr_cutoff
        sig_x[is.na(sig_x)] <- FALSE
    }

    if (!is.null(fdr_y_col)) {
        # FDR filter - if not significant_filter -> FALSE
        sig_y <- df_norm[[fdr_y_col]] <= fdr_cutoff
        sig_y[is.na(sig_y)] <- FALSE
    }
    significant_filter = rep(TRUE, n)
    if (filter_fdr_x) {
        significant_filter = significant_filter & sig_x
    }
    if (filter_fdr_y) {
        significant_filter = significant_filter & sig_y
    }
    sel_x_str <- ifelse(sel_x_cat == "neg",  paste0(x_col, " depleted"),
                   ifelse(sel_x_cat == "pos", paste0(x_col, " enriched"),
                                         paste0(x_col, " no change")))
    sel_y_str <- ifelse(sel_y_cat == "neg",  paste0(y_col, " enriched"),
                   ifelse(sel_y_cat == "pos", paste0(y_col, " enriched"),
                                         paste0(y_col, " no change")))


    if (!is.null(fdr_x_col)) {
        sel_x_str <- ifelse(sig_x,  paste0(sel_x_str, " (sign.)"), paste0(sel_x_str, " (not sign.)"))
    }
    if (!is.null(fdr_y_col)) {
        sel_y_str <- ifelse(sig_y,  paste0(sel_y_str, " (sign.)"), paste0(sel_y_str, " (not sign.)"))
    }
    # verbose hit_class description
    description <- paste(sel_x_str, sel_y_str, diag_desc, sep = " | ")

    # combined hit_class = diagonal + quadrant
    hit_class <- paste(as.character(diagonal), as.character(quadrant), sep = "_")
    hit_class <- factor(hit_class)

    # put all in frame
    df_norm$diagonal    <- diagonal
    df_norm$quadrant    <- quadrant
    df_norm$description <- description
    df_norm$hit_class   <- hit_class
    df_norm$significant_filter <- significant_filter

    # select hits (select + groups)

    # select - filter on diagonal
    if (!is.null(select)) {
        select <- match.arg(select, c("positive", "negative", "both", "none"))
        diag_keep <- switch(
            select,
            "positive" = (diagonal == "pos"),
            "negative" = (diagonal == "neg"),
            "both"     = (diagonal %in% c("pos", "neg")),
            "none"     = (diagonal == "none")
        )
    } else {
        diag_keep <- rep(TRUE, n)  # no filter on diagonal
    }

    # filter on FDR
    if (!is.null(significant_filter)) {
        sig_keep <- significant_filter
    } else {
        sig_keep <- rep(TRUE, n)  # no filter on diagonal
    }

    # group - filter on quadrants (just exact quadrants, no "top"/"mid" etc.)
    if (is.null(groups) || length(groups) == 0 || (!filter_groups)) {
        group_keep <- rep(TRUE, n)
    } else {
        group_keep <- quadrant %in% groups
    }
    selection_idx <- diag_keep & group_keep & sig_keep
    # diese Gene werden gelabelt
    top_labels <- NULL
    # toplabels kann ein logischer Skalar (TRUE/FALSE) oder ein character-Vektor sein.
    if (is.logical(toplabels) && length(toplabels) == 1) {
        if (!is.na(toplabels) && toplabels) {
            top_labels <- as.character(df_norm[[gene_col]][selection_idx])
        } else {
            top_labels <- NULL
        }
    } else if (is.character(toplabels)) {
        top_labels <- toplabels
    } else {
        print(toplabels)
        stop("'toplabels' musst be a length-1 logical or a character vector.")
    }
        
    # ScatterView
    plot_df <- df_norm
    plot_df$X_tmp <- df_norm[[x_col]]
    plot_df$Y_tmp <- df_norm[[y_col]]

    # standard: alle genes are labeled
    plot_df$label_for_plot <- as.character(df_norm[[gene_col]])

    # if label_selected_only = TRUE -> only selected genes are labeled
    if (isTRUE(label_selected_only)) {
        plot_df$label_for_plot[!selection_idx] <- ""
    }

    # if label_selected_only = TRUE, no limit from "top"
    # top_for_plot <- if (isTRUE(label_selected_only)) n else top
    if (is.null(x_cut)){
        x_cut = c(neg_effect_cutoff, pos_effect_cutoff)
    }
    if (is.null(y_cut)){
        y_cut = c(neg_effect_cutoff, pos_effect_cutoff)
    }
    if (data_type == "mle") {
        scatterview <- MAGeCKFlute::ScatterView(
            plot_df,
            x = "X_tmp",
            y = "Y_tmp",
            top = top,
            label = "label_for_plot",
            display_cut = TRUE,
            auto_cut_x = auto_cut_x,
            auto_cut_y = auto_cut_y,
            auto_cut_diag = auto_cut_diag,
            xlab = xlab,
            ylab = ylab,
            groups = groups,
            toplabels = top_labels
        )
    } else {
        scatterview <- MAGeCKFlute::ScatterView(
            plot_df,
            x = "X_tmp",
            y = "Y_tmp",
            model = "custom",
            top = top,
            label = "label_for_plot",
            display_cut = TRUE,
            x_cut = x_cut,
            y_cut = y_cut,
            auto_cut_diag = auto_cut_diag,
            xlab = xlab,
            ylab = ylab,
            slope = slope_fit,
            groups = groups,
            toplabels = top_labels
        )
    }

    # QC fit
    qc_df <- data.frame(
        x = df_norm[[x_col]],
        y = df_norm[[y_col]]
    )

    qc_plot <- ggplot2::ggplot(qc_df, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(alpha = 0.4) +
        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red") +
        ggplot2::labs(
            title = sprintf("Linear fit: %s vs %s", ylab, xlab),
            x = xlab,
            y = ylab
        ) +
        ggplot2::theme_bw()

    # write files
    plot_file      <- file.path(output_dir, paste0(filebase_name, ".jpeg"))
    data_file      <- file.path(output_dir, paste0(filebase_name, "_data.tsv"))
    qc_plot_file   <- file.path(output_dir, paste0(filebase_name, "_lmfit.jpeg"))
    selection_file <- file.path(output_dir, paste0(filebase_name, "_hits_selected.tsv"))

    ggplot2::ggsave(
        filename = plot_file,
        plot = scatterview,
        width = jpeg_width,
        height = jpeg_height
    )

    utils::write.table(
        scatterview$data,
        file = data_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    # write only selected genes, by select + groups
    utils::write.table(
        df_norm[selection_idx, ],
        file = selection_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    ggplot2::ggsave(
        filename = qc_plot_file,
        plot = qc_plot,
        width = jpeg_width,
        height = jpeg_height
    )

    # return this
    res <- list(
        data_type      = data_type,
        scatterview    = scatterview,
        qc_plot        = qc_plot,
        fit            = fit,
        slope_fit      = slope_fit,
        delta_cutoff   = delta_cutoff,
        intercept_min  = intercept_min,
        intercept_max  = intercept_max,
        data           = df_norm,
        plot_file      = plot_file,
        data_file      = data_file,
        qc_plot_file   = qc_plot_file,
        hits_file      = selection_file
    )

    invisible(res)
}

# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
# When called as  Rscript mageck.R <params.json>  the block below runs;
# when the file is sourced interactively it is skipped.
if (!interactive()) {
  suppressPackageStartupMessages(library(jsonlite))

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) {
    stop("Usage: Rscript mageck.R <params.json>", call. = FALSE)
  }

  p <- fromJSON(args[1], simplifyVector = FALSE)

  # Mandatory fields
  mode <- p$mode %||% "scatter_view"

  if (mode == "scatter_view") {
    res <- RunMageckScatterView(
      input_file         = p$input_file,
      data_type          = p$data_type          %||% "mle",
      output_dir         = p$output_dir,
      filebase_name      = p$filebase_name,
      x_col              = p$x_col,
      y_col              = p$y_col,
      gene_col           = p$gene_col           %||% "id",
      sep                = p$sep                %||% "\t",
      normalize          = p$normalize          %||% TRUE,
      top                = p$top                %||% 10L,
      groups             = unlist(p$groups)     %||% c("bottomleft"),
      select             = p$select,
      neg_effect_cutoff  = p$neg_effect_cutoff  %||% -0.4,
      pos_effect_cutoff  = p$pos_effect_cutoff  %||%  0.4,
      delta_cutoff_k     = p$delta_cutoff_k     %||%  2,
      filter_fdr_x       = p$filter_fdr_x       %||% FALSE,
      filter_fdr_y       = p$filter_fdr_y       %||% FALSE,
      filter_groups      = p$filter_groups      %||% TRUE,
      fdr_cutoff         = p$fdr_cutoff         %||% 0.05,
      toplabels          = p$toplabels          %||% TRUE,
      label_selected_only = p$label_selected_only %||% FALSE,
      xlab               = p$xlab,
      ylab               = p$ylab,
      jpeg_width         = p$jpeg_width         %||% 20,
      jpeg_height        = p$jpeg_height        %||% 15,
      auto_cut_diag      = p$auto_cut_diag      %||% 2,
      auto_cut_x         = p$auto_cut_x         %||% 2,
      auto_cut_y         = p$auto_cut_y         %||% 2,
      x_cut              = unlist(p$x_cut),
      y_cut              = unlist(p$y_cut)
    )
    message("MAGeCK scatter_view done.")
    message("  plot     -> ", res$plot_file)
    message("  data     -> ", res$data_file)
    message("  qc_plot  -> ", res$qc_plot_file)
    message("  hits     -> ", res$hits_file)
  } else {
    stop(paste("Unknown mode:", mode), call. = FALSE)
  }
}
