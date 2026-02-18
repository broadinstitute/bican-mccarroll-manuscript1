# library (ggplot2)
# library (cowplot)

age_preds_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/age_prediction/age_prediction_region_alpha_0/age_prediction_results_alpha0_donor_predictions.txt"
cell_type_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/metadata/cell_types_for_prs_test.txt"

# summarize median absolute error, mean absolute error, and median absolute error
# in the youngest 20% of donors for each cell type and region, and save to a summary file

#' Summarize age-prediction model performance and covariate associations
#'
#' This function reads per-donor age prediction results, optionally restricts
#' to a user-specified set of cell types, computes model-level error metrics
#' for each cell type by region, generates summary plots, and fits simple
#' regression models relating prediction error to model size / data size
#' covariates.
#'
#' The expected input file is a tab-delimited table with one row per donor per
#' cell type x region model. At minimum, it should contain the columns needed
#' by \code{\link{compute_age_error_metrics}} (including \code{cell_type},
#' \code{region}, \code{donor}, \code{age}, the prediction column used by
#' \code{compute_age_error_metrics}, and \code{num_nuclei}, \code{num_umis},
#' \code{num_features}).
#'
#' @param age_preds_file Character scalar. Path to a tab-delimited file
#' containing per-donor predictions and metadata.
#' @param cell_type_file Optional character scalar. Path to a text file
#' containing one cell type per line. If provided, \code{age_preds_file} is
#' filtered to those cell types.
#' @param out_summary_file Character scalar. Path for writing a tab-delimited
#' summary table of model-level metrics (one row per cell type x region).
#' @param out_pdf_file Optional character scalar. If provided, summary plots
#' are written to this PDF.
#'
#' @return The list returned by \code{\link{summarize_age_prediction_results}}.
#'
#' @export
summarize_age_prediction_results_file <- function(age_preds_file,
                                                  cell_type_file = NULL,
                                                  out_summary_file,
                                                  out_pdf_file = NULL) {
    age_preds <- utils::read.table(
        age_preds_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )

    cell_types <- NULL
    if (!is.null(cell_type_file)) {
        cell_types <- utils::read.table(
            cell_type_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE
        )$V1
        cell_types <- as.character(cell_types)
    }

    res <- summarize_age_prediction_results(
        age_preds = age_preds,
        cell_types = cell_types
    )

    utils::write.table(
        res$metrics_df, file = out_summary_file, sep = "\t",
        quote = FALSE, row.names = FALSE, col.names = TRUE
    )

    if (!is.null(out_pdf_file)) {
        y_axis_label_mae <- cowplot::ggdraw() +
            cowplot::draw_label("Mean AE (decades)", angle = 90, size = 13)

        combined <- cowplot::plot_grid(
            y_axis_label_mae,
            cowplot::plot_grid(
                res$plots$p_feat_mae,
                res$plots$p_nuc_mae,
                res$plots$p_umi_mae,
                ncol = 1
            ),
            ncol = 2, rel_widths = c(0.05, 1)
        )

        y_axis_label_y20 <- cowplot::ggdraw() +
            cowplot::draw_label("Mean AE in youngest 20% (decades)", angle = 90, size = 13)

        combined_y20 <- cowplot::plot_grid(
            y_axis_label_y20,
            cowplot::plot_grid(
                res$plots$p_feat_y20,
                res$plots$p_nuc_y20,
                res$plots$p_umi_y20,
                ncol = 1
            ),
            ncol = 2, rel_widths = c(0.05, 1)
        )

        grDevices::pdf(out_pdf_file, width = 10, height = 7)
        on.exit(grDevices::dev.off(), add = TRUE)

        print(res$plots$p_mae)
        print(res$plots$p_mae_y20)
        print(combined)
        print(combined_y20)
    }

    res
}

#' Summarize age-prediction model performance and covariate associations
#'
#' This function analyzes per-donor age prediction results without performing
#' any file I/O. It optionally restricts to a user-specified set of cell types,
#' computes model-level error metrics for each cell type by region, generates
#' summary plots, and fits simple regression models relating prediction error
#' to model size / data size covariates.
#'
#' @param age_preds data.frame containing per-donor predictions and metadata.
#' @param cell_types Optional character vector of cell types to keep. If
#' provided, \code{age_preds} is filtered to those cell types.
#'
#' @details
#' The returned \code{plots} element contains individual plots:
#' \describe{
#'   \item{p_mae}{Heatmap of \code{mae}.}
#'   \item{p_mae_y20}{Heatmap of \code{mae_young20}.}
#'   \item{p_feat_mae}{MAE vs \code{num_features}.}
#'   \item{p_nuc_mae}{MAE vs \code{num_nuclei}.}
#'   \item{p_umi_mae}{MAE vs \code{num_umis}.}
#'   \item{p_feat_y20}{MAE in youngest 20\% vs \code{num_features}.}
#'   \item{p_nuc_y20}{MAE in youngest 20\% vs \code{num_nuclei}.}
#'   \item{p_umi_y20}{MAE in youngest 20\% vs \code{num_umis}.}
#' }
#'
#' Combined panels are intended to be constructed by the file I/O wrapper.
#'
#' @export
summarize_age_prediction_results <- function(age_preds, cell_types = NULL) {
    stopifnot(is.data.frame(age_preds))

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }
    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("Package 'cowplot' is required.")
    }

    if (!is.null(cell_types)) {
        if (!is.character(cell_types)) {
            stop("'cell_types' must be a character vector.")
        }
        age_preds <- age_preds[age_preds$cell_type %in% cell_types, , drop = FALSE]
    }

    metrics_df <- compute_age_error_metrics(age_preds)

    p_mae <- plot_age_metric_heatmap(metrics_df, metric = "mae", text_size = 4)
    p_mae_y20 <- plot_age_metric_heatmap(metrics_df, metric = "mae_young20")

    fit_all <- fit_age_mae_regression(metrics_df, outcome = "mae")
    fit_y20 <- fit_age_mae_regression(metrics_df, outcome = "mae_young20")

    p_feat_mae <- plot_error_vs_predictor(
        metrics_df, x_col = "num_features", y_col = "mae",
        x_scale = "identity", annotate_r2_size = 5, show_y_label = FALSE
    ) +  ggplot2::labs(x = "Age-associated DE genes")

    p_nuc_mae <- plot_error_vs_predictor(
        metrics_df, x_col = "num_nuclei", y_col = "mae",
        x_scale = "identity", annotate_r2_size = 5, show_y_label = FALSE
    ) +  ggplot2::labs(x = "Total Nuclei")

    p_umi_mae <- plot_error_vs_predictor(
        metrics_df, x_col = "num_umis", y_col = "mae",
        x_scale = "identity", annotate_r2_size = 5, show_y_label = FALSE
    ) + ggplot2::labs(x = "Total UMIs")

    p_feat_y20 <- plot_error_vs_predictor(
        metrics_df, x_col = "num_features", y_col = "mae_young20",
        x_scale = "identity", annotate_r2_size = 5, show_y_label = FALSE
    ) +  ggplot2::labs(x = "Age-associated DE genes")

    p_nuc_y20 <- plot_error_vs_predictor(
        metrics_df, x_col = "num_nuclei", y_col = "mae_young20",
        x_scale = "identity", annotate_r2_size = 5, show_y_label = FALSE
    ) +  ggplot2::labs(x = "Total Nuclei")

    p_umi_y20 <- plot_error_vs_predictor(
        metrics_df, x_col = "num_umis", y_col = "mae_young20",
        x_scale = "identity", annotate_r2_size = 5, show_y_label = FALSE
    ) + ggplot2::labs(x = "Total UMIs")


    plots <- list(
        p_mae = p_mae,
        p_mae_y20 = p_mae_y20,
        p_feat_mae = p_feat_mae,
        p_nuc_mae = p_nuc_mae,
        p_umi_mae = p_umi_mae,
        p_feat_y20 = p_feat_y20,
        p_nuc_y20 = p_nuc_y20,
        p_umi_y20 = p_umi_y20
    )

    list(
        age_preds = age_preds,
        metrics_df = metrics_df,
        fit_all = fit_all,
        fit_y20 = fit_y20,
        plots = plots
    )
}







#' Summarize age-prediction error metrics per cell_type x region
#'
#' @param age_preds data.frame with per-donor predictions and metadata.
#' Required columns: cell_type, region, donor, age, pred_mean_corrected,
#' num_nuclei, num_umis, num_features.
#' @param group_cols character vector of grouping columns.
#' @param pred_col prediction column (corrected predictions).
#' @param age_col chronological age column.
#' @param young_q quantile defining "youngest" subset within each group.
#'
#' @return data.frame with one row per group and error metrics/covariates.
compute_age_error_metrics <- function(age_preds,
                                      group_cols = c("cell_type", "region"),
                                      pred_col = "pred_mean",
                                      age_col = "age",
                                      young_q = 0.2) {
    stopifnot(is.data.frame(age_preds))
    req <- c(group_cols, "donor", age_col, pred_col,
             "num_nuclei", "num_umis", "num_features")
    miss <- setdiff(req, colnames(age_preds))
    if (length(miss) > 0) {
        stop("Missing required columns: ", paste(miss, collapse = ", "))
    }

    dt <- data.table::as.data.table(age_preds)

    #Make R CMD CHECK Happy
    abs_err <- young_cut <- is_young <- num_features <- NULL

    dt[, abs_err := abs(get(pred_col) - get(age_col))]

    dt[, young_cut := stats::quantile(get(age_col), probs = young_q, type = 7, na.rm = TRUE),
       by = group_cols]
    dt[, is_young := get(age_col) <= young_cut]

    summarize_one_group <- function(d) {
        nf_u <- unique(d[!is.na(num_features), num_features])
        nf <- if (length(nf_u) == 1) as.numeric(nf_u) else NA_real_

        if (length(nf_u) > 1) {
            grp_vals <- vapply(group_cols, function(gc) as.character(d[[gc]][1]), character(1))
            warning("num_features not unique within group: ",
                    paste(grp_vals, collapse = " / "))
        }

        # Enforce numeric/double outputs where appropriate to keep data.table happy
        list(
            n_donors = as.integer(data.table::uniqueN(d$donor)),
            mae = as.numeric(mean(d$abs_err, na.rm = TRUE)),
            medae = as.numeric(stats::median(d$abs_err, na.rm = TRUE)),
            mae_young20 = as.numeric(mean(d$abs_err[d$is_young], na.rm = TRUE)),
            medae_young20 = as.numeric(stats::median(d$abs_err[d$is_young], na.rm = TRUE)),
            num_nuclei = as.numeric(sum(d$num_nuclei, na.rm = TRUE)),
            num_umis = as.numeric(sum(d$num_umis, na.rm = TRUE)),
            num_features = as.numeric(nf)
        )
    }

    out <- dt[, summarize_one_group(.SD), by = group_cols]
    data.table::setDF(out)
}

#' Compute pairwise residual age correlations within and across cell populations
#'
#' Computes pairwise correlations of corrected residual age values
#' (`resid_mean_corrected`) across two comparison schemes:
#' (1) between cell types within each region and
#' (2) between regions within each cell type.
#' The function returns two numeric vectors containing all pairwise
#' correlation coefficients for each comparison class.
#'
#' @param model_predictions A data frame containing corrected residual
#' age values (`resid_mean_corrected`) along with `cell_type`, `region`,
#' and `donor` columns.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{within_region}{Pairwise correlations across cell types within regions.}
#'   \item{within_celltype}{Pairwise correlations across regions within cell types.}
#' }
#'
#' @export
compute_residual_age_correlations <- function(model_predictions) {

    df <- model_predictions

    # Keep only required columns
    df <- df[, c("cell_type", "region", "donor", "resid_mean_corrected")]

    # ----------------------------
    # 1) Within-region, across cell types
    # ----------------------------
    regions <- unique(df$region)
    cell_types <- unique(df$cell_type)

    within_region <- c()

    for (r in regions) {

        df_r <- df[df$region == r, ]

        ct_levels <- unique(df_r$cell_type)

        if (length(ct_levels) < 2) next

        combs <- utils::combn(ct_levels, 2, simplify = FALSE)

        for (pair in combs) {

            df1 <- df_r[df_r$cell_type == pair[1], ]
            df2 <- df_r[df_r$cell_type == pair[2], ]

            merged <- merge(df1[, c("donor", "resid_mean_corrected")],
                            df2[, c("donor", "resid_mean_corrected")],
                            by = "donor",
                            suffixes = c("_1", "_2"))

            if (nrow(merged) >= 10) {
                cor_val <- cor(merged$resid_mean_corrected_1,
                               merged$resid_mean_corrected_2,
                               use = "complete.obs")
                within_region <- c(within_region, cor_val)
            }
        }
    }

    # ----------------------------
    # 2) Within-cell-type, across regions
    # ----------------------------
    within_celltype <- c()

    for (ct in cell_types) {

        df_ct <- df[df$cell_type == ct, ]

        region_levels <- unique(df_ct$region)

        if (length(region_levels) < 2) next

        combs <- utils::combn(region_levels, 2, simplify = FALSE)

        for (pair in combs) {

            df1 <- df_ct[df_ct$region == pair[1], ]
            df2 <- df_ct[df_ct$region == pair[2], ]

            merged <- merge(df1[, c("donor", "resid_mean_corrected")],
                            df2[, c("donor", "resid_mean_corrected")],
                            by = "donor",
                            suffixes = c("_1", "_2"))

            if (nrow(merged) >= 10) {
                cor_val <- cor(merged$resid_mean_corrected_1,
                               merged$resid_mean_corrected_2,
                               use = "complete.obs")
                within_celltype <- c(within_celltype, cor_val)
            }
        }
    }

    return(list(
        within_region = within_region,
        within_celltype = within_celltype
    ))
}

#' Plot distributions of residual age correlations
#'
#' Generates a violin plot with embedded boxplots comparing the
#' distribution of residual age correlations computed within regions
#' (across cell types) and within cell types (across regions).
#' This visualization is intended to assess whether donor-specific
#' variation in residual age is more consistent within a cell type
#' across regions than across cell types within a region.
#'
#' @param corr_results A list returned by
#' `compute_residual_age_correlations()`, containing numeric vectors
#' `within_region` and `within_celltype`.
#' @param fill_colors A named character vector specifying fill colors
#' for the two comparison groups.
#'
#' @return A ggplot object.
#'
#' @export
plot_residual_age_correlation_distributions <- function(
        corr_results,
        fill_colors = c(
            "Within region\n(across cell types)" = "#8DA0CB",
            "Within cell type\n(across regions)" = "#66C2A5"
        )
) {

    correlation <- group <- NULL

    df_plot <- data.frame(
        correlation = c(corr_results$within_region,
                        corr_results$within_celltype),
        group = factor(
            c(rep("Within region\n(across cell types)",
                  length(corr_results$within_region)),
              rep("Within cell type\n(across regions)",
                  length(corr_results$within_celltype))),
            levels = c("Within region\n(across cell types)",
                       "Within cell type\n(across regions)")
        )
    )

    ggplot2::ggplot(df_plot,
                    ggplot2::aes(x = group, y = correlation, fill = group)
    ) +
        ggplot2::geom_violin(alpha = 0.6, width = 0.8, color = NA) +
        ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA,
                              alpha = 0.8, color = "black") +
        ggplot2::scale_fill_manual(values = fill_colors) +
        ggplot2::labs(
            x = NULL,
            y = "Residual age correlation"
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "none")
}




#' Heatmap for a single metric across cell_type x region
#'
#' Creates a tile heatmap of one column from `metrics_df` (e.g. MAE) across
#' cell types (rows) and regions (columns). Optionally overlays formatted
#' values as text, with text color chosen for readability against the fill.
#'
#' @param metrics_df Data.frame returned by `compute_age_error_metrics()`.
#' @param metric Character scalar giving the column name to plot (e.g. `"mae"`
#'   or `"mae_young20"`).
#' @param cell_type_col Character scalar giving the column name in `metrics_df`
#'   that contains cell type.
#' @param region_col Character scalar giving the column name in `metrics_df`
#'   that contains region.
#' @param na_color Color used for missing tiles (passed as `na.value` to the
#'   fill scale).
#' @param add_text Logical; if `TRUE`, overlay numeric values on tiles.
#' @param text_digits Integer number of digits after the decimal to display
#'   when `add_text = TRUE`.
#' @param text_size Numeric text size for the overlay labels when
#'   `add_text = TRUE`.
#' @param text_threshold Optional numeric threshold used to choose label text
#'   color. Values less than or equal to this threshold use white text and
#'   values greater than this threshold use black text. If `NULL`, a threshold
#'   is computed from the data range.
#'
#' @return A ggplot object.
#' @export
plot_age_metric_heatmap <- function(metrics_df,
                                metric,
                                cell_type_col = "cell_type",
                                region_col = "region",
                                na_color = "grey90",
                                add_text = TRUE,
                                text_digits = 1,
                                text_size = 3,
                                text_threshold = NULL) {
    stopifnot(is.data.frame(metrics_df))
    need <- c(cell_type_col, region_col, metric)
    miss <- setdiff(need, colnames(metrics_df))
    if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

    #Make R CMD CHECK Happy
    cell_type <- region <- value <- label <- text_col <- cell_type_label <- region_label <- mean_value <- NULL

    dt <- data.table::as.data.table(metrics_df)
    dt[, value := get(metric)]
    data.table::setnames(dt, c(cell_type_col, region_col), c("cell_type", "region"))

    # Pretty labels (for axes / legend)
    dt[, cell_type_label := gsub("_", " ", cell_type, fixed = TRUE)]
    dt[, region_label := gsub("_", " ", region, fixed = TRUE)]
    metric_label <- gsub("_", " ", metric, fixed = TRUE)

    # Sort rows by average metric across regions (best = lowest error at top)
    #cell_means <- dt[, .(mean_value = mean(value, na.rm = TRUE)), by = cell_type_label]
    cell_means <- dt[, list(mean_value = mean(value, na.rm = TRUE)), by = cell_type_label]
    cell_means <- cell_means[!is.nan(mean_value)]
    # For errors, lower is better; to put best at top, make it the last factor level
    cell_levels <- cell_means[order(mean_value, decreasing = TRUE), cell_type_label]
    dt[, cell_type_label := factor(cell_type_label, levels = cell_levels)]

    # Choose text color based on fill intensity:
    # With "Blues 3" here, low values are darker -> use white text; high values are lighter -> black text.
    rng <- range(dt$value, na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
        # degenerate case
        dt[, text_col := "black"]
    } else {
        if (is.null(text_threshold)) {
            text_threshold <- rng[1] + 0.55 * (rng[2] - rng[1])
        }
        dt[, text_col := ifelse(value <= text_threshold, "white", "black")]
    }

    p <- ggplot2::ggplot(dt, ggplot2::aes(x = region_label, y = cell_type_label, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradientn(
            colors = grDevices::hcl.colors(256, "Blues 3"),
            na.value = na_color
        ) +
        ggplot2::labs(x = "Region", y = "Cell type", fill = metric_label) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank()
        )

    if (isTRUE(add_text)) {
        dt_txt <- dt[!is.na(value)]
        dt_txt[, label := sprintf(paste0("%.", text_digits, "f"), value)]
        p <- p + ggplot2::geom_text(
            data = dt_txt,
            ggplot2::aes(label = label, color = text_col),
            size = text_size,
            show.legend = FALSE
        ) +
            ggplot2::scale_color_identity()
    }

    p
}


fit_age_mae_regression <- function(metrics_df,
                                            outcome,
                                            numeric_covars = c("num_nuclei",
                                                               "num_umis",
                                                               "num_features"),
                                            log_pseudocount = 1) {
    stopifnot(is.data.frame(metrics_df))

    need <- c(outcome, numeric_covars)
    miss <- setdiff(need, colnames(metrics_df))
    if (length(miss) > 0) {
        stop("Missing required columns: ", paste(miss, collapse = ", "))
    }

    df <- metrics_df

    # Transform numeric predictors: scale(log10(x))
    for (nm in numeric_covars) {
        x <- df[[nm]]
        x2 <- log10(pmax(x, log_pseudocount, na.rm = TRUE))
        df[[paste0("zlog10_", nm)]] <- as.numeric(scale(x2))
    }

    rhs <- paste(paste0("zlog10_", numeric_covars), collapse = " + ")
    fml <- stats::as.formula(paste(outcome, "~", rhs))

    fit <- stats::lm(fml, data = df)

    sm <- summary(fit)
    coefs <- sm$coefficients
    ci <- try(stats::confint(fit), silent = TRUE)

    coef_df <- data.frame(
        term = rownames(coefs),
        estimate = coefs[, "Estimate"],
        std.error = coefs[, "Std. Error"],
        statistic = coefs[, "t value"],
        p.value = coefs[, "Pr(>|t|)"],
        stringsAsFactors = FALSE
    )

    if (!inherits(ci, "try-error")) {
        coef_df$conf.low <- ci[coef_df$term, 1]
        coef_df$conf.high <- ci[coef_df$term, 2]
    } else {
        coef_df$conf.low <- NA_real_
        coef_df$conf.high <- NA_real_
    }

    list(
        fit = fit,
        coef_table = coef_df,
        r.squared = sm$r.squared,
        adj.r.squared = sm$adj.r.squared
    )
}




plot_mae_vs_num_features <- function(metrics_df,
                                     x_col = "num_features",
                                     y_col = "mae",
                                     color_by = NULL,
                                     smoother_span = 1,
                                     point_size = 2.6,
                                     alpha = 0.9,
                                     show_legend = TRUE) {

    stopifnot(is.data.frame(metrics_df))

    need <- c(x_col, y_col)
    if (!is.null(color_by)) need <- c(need, color_by)

    miss <- setdiff(need, colnames(metrics_df))
    if (length(miss) > 0) {
        stop("Missing required columns: ", paste(miss, collapse = ", "))
    }

    num_features <- mae <- grp <- NULL

    df <- metrics_df
    df$num_features <- df[[x_col]]
    df$mae <- df[[y_col]]

    if (!is.null(color_by)) {
        raw_grp <- as.factor(df[[color_by]])
        lbls <- gsub("_", " ", levels(raw_grp))
        df$grp <- factor(raw_grp, levels = levels(raw_grp), labels = lbls)
    }

    n_grp <- length(levels(df$grp))
    pal <- .make_palette(n_grp)

    if (is.null(color_by)) {
        p <- ggplot2::ggplot(df, ggplot2::aes(x = num_features, y = mae)) +
            ggplot2::geom_point(
                size = point_size,
                alpha = alpha,
                shape = 21,
                fill = "grey70",
                color = "black",
                stroke = 0.6
            )
    } else {
        n_grp <- length(levels(df$grp))
        if (n_grp > length(pal)) {
            stop("Too many groups in ", color_by, " (", n_grp, "). Increase palette length.")
        }

        p <- ggplot2::ggplot(df, ggplot2::aes(x = num_features, y = mae, fill = grp)) +
            ggplot2::geom_point(
                size = point_size,
                alpha = alpha,
                shape = 21,
                color = "black",
                stroke = 0.5
            ) +
            ggplot2::scale_fill_manual(values = pal[seq_len(n_grp)], name = color_by)
    }

    p <- p +
        ggplot2::geom_smooth(
            method = "loess",
            span = smoother_span,
            se = FALSE,
            color = "grey25",
            linewidth = 1.1,
            inherit.aes = FALSE,
            mapping = ggplot2::aes(x = num_features, y = mae),
            data = df
        ) +
        ggplot2::labs(
            x = "Number of age-informative genes",
            y = "Mean absolute error (years)"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            axis.title = ggplot2::element_text(size = 13),
            axis.text = ggplot2::element_text(size = 11),
            legend.position = if (show_legend && !is.null(color_by)) "right" else "none",
            legend.title = ggplot2::element_text(size = 12),
            legend.text = ggplot2::element_text(size = 10)
        )

    p
}

.make_palette <- function(n) {
    if (requireNamespace("colorspace", quietly = TRUE)) {
        # "Dark 3" can collide; "Set 3" often includes pale/yellowish tones.
        # "Glasbey" is ideal but not always available. This is a good compromise.
        colorspace::qualitative_hcl(n, palette = "Dark 2")
    } else {
        # Hand-curated, distinct, no yellow, no black/grey.
        c(
            "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
            "#66A61E", "#E6AB02", "#A6761D", "#666666",
            "#1F78B4", "#B2DF8A", "#FB9A99", "#CAB2D6"
        )[seq_len(n)]
    }
}

plot_error_vs_predictor <- function(metrics_df,
                                    x_col,
                                    y_col = "mae",
                                    x_label = NULL,
                                    y_label = NULL,
                                    x_scale = c("identity", "log10"),
                                    smoother_span = 1,
                                    add_smoother = TRUE,
                                    annotate_adj_r2 = TRUE,
                                    point_size = 2.6,
                                    alpha = 0.9,
                                    annotate_r2_size = 3.6,
                                    log_pseudocount = 1,
                                    show_y_label = TRUE) {

    stopifnot(is.data.frame(metrics_df))
    x_scale <- match.arg(x_scale)

    need <- c(x_col, y_col)
    miss <- setdiff(need, colnames(metrics_df))
    if (length(miss) > 0) {
        stop("Missing required columns: ", paste(miss, collapse = ", "))
    }

    num_x <- y_val <- NULL

    df <- metrics_df
    df$num_x <- df[[x_col]]
    df$y_val <- df[[y_col]]

    if (is.null(x_label)) {
        x_label <- gsub("_", " ", x_col)
    } else {
        x_label <- gsub("_", " ", x_label)
    }

    if (is.null(y_label)) {
        if (tolower(y_col) == "mae") {
            y_label <- "Mean AE (years)"
        } else if (tolower(y_col) == "medae") {
            y_label <- "Median AE (years)"
        } else if (tolower(y_col) == "mae_young20") {
            y_label <- "Mean AE in youngest 20% (years)"
        } else if (tolower(y_col) == "medae_young20") {
            y_label <- "Median AE in youngest 20% (years)"
        } else {
            y_label <- gsub("_", " ", y_col)
        }
    } else {
        y_label <- gsub("_", " ", y_label)
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = num_x, y = y_val)) +
        ggplot2::geom_point(
            size = point_size,
            alpha = alpha,
            shape = 21,
            fill = "grey70",
            color = "black",
            stroke = 0.5
        ) +
        ggplot2::labs(
            x = x_label,
            y = if (isTRUE(show_y_label)) y_label else NULL
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            axis.title = ggplot2::element_text(size = 13),
            axis.text = ggplot2::element_text(size = 11)
        )

    if (x_scale == "log10") {
        p <- p + ggplot2::scale_x_log10()
    }

    if (add_smoother) {
        p <- p + ggplot2::geom_smooth(
            method = "loess",
            span = smoother_span,
            se = FALSE,
            color = "grey25",
            linewidth = 1.1
        )
    }

    if (annotate_adj_r2) {
        x_raw <- df$num_x
        x2 <- log10(pmax(x_raw, log_pseudocount, na.rm = TRUE))
        z_x <- as.numeric(scale(x2))

        fit_df <- data.frame(y = df$y_val, z_x = z_x)
        fit <- stats::lm(y ~ z_x, data = fit_df)
        adj_r2 <- summary(fit)$adj.r.squared

        xr <- range(df$num_x, na.rm = TRUE)
        yr <- range(df$y_val, na.rm = TRUE)
        ann_x <- xr[1] + 0.70 * (xr[2] - xr[1])
        ann_y <- yr[1] + 0.92 * (yr[2] - yr[1])

        p <- p + ggplot2::annotate(
            "text",
            x = ann_x,
            y = ann_y,
            label = sprintf("adj R2 = %.2f", adj_r2),
            hjust = 0,
            size = annotate_r2_size,
            color = "grey20"
        )
    }

    p
}
